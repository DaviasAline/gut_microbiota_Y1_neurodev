## Aline Davias
## 21/02/2024
## Analyses de la beta diversité 

# Packages, functions and data loading ----
library(tidyverse)
library(vegan)
library(phyloseq)
library(broom)
library(forestplot)
library(bkmr)
library(fields)
library(psych)
library(expss)
library(SRS)
library(see)
library(writexl)
library(forcats)
library(egg)
library(corrplot)
load("1_intermediate_data/3_data_imputation_AD_gumme.RData")
source("3_programs/4_functions_AD_gumme.R")
rm(bdd_final_imp_1_sensi_seuil, bdd_final_imp_sensi_seuil, bdd_genera_imp, 
   covar_a_imputer, covar_a_tester, 
   comp_effectifs, model, model_sensi_home, model_sensi_seuil, 
   tbl_model_covar, tbl_model_final, 
   table_cor, table_cor_sg, model_Y3_1, model_Y3_2, heatmap_cor, heatmap_cor_pairwise)
asv_raw_not_rarefied <- read_labelled_csv("0_source_data/gut_microbiota_ASVbased_Y1_labelled_AD_20220504_7239.csv") 

# Vectors ----
covariates <- covar_vec_model_final
outcomes <- bdd %>%
  select(all_of(neuro_vec)) %>%
  select(-ch_socawar_y3,
         -ch_soccog_y3,
         -ch_soccom_y3, 
         -ch_socmot_y3, 
         -ch_RRB_y3) %>%
  colnames()
rm(covar_vec_model_final, neuro_vec, genera_linear, microbiote_vec)

spanner_names <- c("Internalizing CBCL score at 2 years",
                   "Externalizing CBCL score at 2 years",
                   "Total SRS score at 3 years",
                   "Inhibition BRIEF-P score at 3 years",
                   "Shift BRIEF-P score at 3 years",
                   "Emotional control BRIEF-P score at 3 years",
                   "Working memory BRIEF-P score at 3 years",
                   "Plan and organization BRIEF-P score at 3 years",
                   "Verbal comprehension IQ score at 3 years",
                   "Visuospatial IQ score at 3 years",
                   "Work memory IQ score at 3 years",
                   "Total IQ score at 3 years")

# Création variables neuro catégorisées (tertiles) ----
bdd_final_imp_1 <- bdd_final_imp_1 %>%
  mutate(across(all_of(outcomes), 
                ~cut(., 
                     breaks = quantile(., probs = 0:3/3, na.rm = TRUE), 
                     labels = c("1st tertile", "2nd tertile", "3rd tertile"), 
                     include.lowest = TRUE),
                .names = "{.col}_cat"))
outcomes <- paste0(outcomes, "_cat")

bdd_final_imp_1 %>%
  select(all_of(outcomes)) %>%
  tbl_summary()

# Beta diversity ----
## Data cleaning ----
### Rarefaction ----
#### Data preparation 
asv_raw_not_rarefied <- 
  asv_raw_not_rarefied %>%
  select(ident, 
         starts_with("ch_feces_raw_asv")) %>%
  na.omit()
row.names(asv_raw_not_rarefied) <- NULL 
asv_raw_not_rarefied <- asv_raw_not_rarefied %>%
  column_to_rownames("ident") %>%
  t() %>%
  as.data.frame()%>% 
  rownames_to_column("ch_feces_ASV_ID_Y1") %>%
  mutate(
    ch_feces_ASV_ID_Y1 = str_replace_all(ch_feces_ASV_ID_Y1, c("ch_feces_raw_" = "", "_Y1"=""))) %>%
  column_to_rownames("ch_feces_ASV_ID_Y1") %>%
  otu_table(taxa_are_rows = TRUE)     

#### from the raw ASV table, choose a threshold for the sequencing depth 
#### define a subset of samples to keep = samples with a sequencing depth > the chosen threshold  
#### samples with a sequencing depth < the chosen threshold become missing data 
keep_5000 <- 
  names(which(sample_sums(asv_raw_not_rarefied)>= 5000)) %>%
  prune_samples(asv_raw_not_rarefied) %>% 
  as.data.frame()                # Loss of 6 samples 

#### reduce the sequencing depth of the samples to the chosen threshold
#### the sequences kept within each sample are randomly selected 
ASV_rarefied_5000_Y1 <- keep_5000 %>%
  SRS(5000, set_seed = TRUE, seed = 1)
rownames(ASV_rarefied_5000_Y1)<- rownames(keep_5000)

#### put the dataframe with rarefied ASVs in columns and samples in rows
ASV_rarefied_5000_Y1 <- 
  ASV_rarefied_5000_Y1 %>%
  t() %>%
  as.data.frame()
rm(keep_5000)

#### check if the rarefaction worked properly 
#### ok the rowsums are equal to 5000, the threshold we chose
rowSums(ASV_rarefied_5000_Y1)   

### Metric calculation (Bray Curtis) ----
### Calcul de la metric Bray Curtis
set.seed(1996)
metric_bray_curtis <- vegan::vegdist(ASV_rarefied_5000_Y1, method = "bray")
### warning message: Plus d’une classe "dist" est trouvée en cache : Utilisation de la première, depuis l’espace de noms 'BiocGenerics'. Aussi défini par ‘spam’
### ok, phyloseq package is supposed to use BiocGeneric 
nmds <- metaMDS(metric_bray_curtis)
scores(nmds) %>%
  as_tibble(rownames = "ident") %>%
  ggplot(aes(x=NMDS1, y=NMDS2)) +
  geom_point()
# nmds <- metaMDS(ASV_rarefied_5000_Y1, autotransform = FALSE) # ?
# scores(nmds) %>%
#   as_tibble(rownames = "ident") %>%
#   ggplot(aes(x=NMDS1, y=NMDS2)) +
#   geom_point()

metric_bray_curtis <- 
  metric_bray_curtis %>% 
  as.matrix %>% 
  as.data.frame() %>%
  rownames_to_column(var = "ident")

### Metadata ----
metadata <- bdd_final_imp_1 %>%
  select(ident, 
         all_of(covariates),
         all_of(covar_age),
         all_of(outcomes)) %>%
  mutate(ident = as.character(ident))

#### merge betadiversity data and metadata
bdd_final <- inner_join(metric_bray_curtis, metadata, by = "ident")

#### filter CBCL_int (because of the NA) ----
excluded_cbcl_int <- bdd_final %>% select(ident, ch_cbclintscore_y2_cat) %>% filter_all(any_vars(is.na(.)))
excluded_cbcl_int <- excluded_cbcl_int$ident

bdd_final_cbcl_int <- bdd_final %>% 
  select(-all_of(excluded_cbcl_int)) %>%
  filter(!ident %in% excluded_cbcl_int)
all_dist_cbcl_int <- bdd_final_cbcl_int %>%
  select(all_of(.[["ident"]])) %>%
  as.dist()

#### filter CBCL ext (because of the NA) ----
excluded_cbcl_ext <- bdd_final %>% select(ident, ch_cbclextscore_y2_cat) %>% filter_all(any_vars(is.na(.)))
excluded_cbcl_ext <- excluded_cbcl_ext$ident

bdd_final_cbcl_ext <- bdd_final %>% 
  select(-all_of(excluded_cbcl_ext)) %>%
  filter(!ident %in% excluded_cbcl_ext)
all_dist_cbcl_ext <- bdd_final_cbcl_ext %>%
  select(all_of(.[["ident"]])) %>%
  as.dist()

#### filter SRS total (because of the NA) ----
excluded_srs <- bdd_final %>% select(ident, ch_SRStotal_y3_cat) %>% filter_all(any_vars(is.na(.)))
excluded_srs <- excluded_srs$ident

bdd_final_srs <- bdd_final %>% 
  select(-all_of(excluded_srs)) %>%
  filter(!ident %in% excluded_srs)
all_dist_srs <- bdd_final_srs %>%
  select(all_of(.[["ident"]])) %>%
  as.dist()

#### filter brief p inhibit (because of the NA) ----
excluded_briefp_inhibit <- bdd_final %>% select(ident, ch_briefpinhibit_y3_cat) %>% filter_all(any_vars(is.na(.)))
excluded_briefp_inhibit <- excluded_briefp_inhibit$ident

bdd_final_briefp_inhibit <- bdd_final %>% 
  select(-all_of(excluded_briefp_inhibit)) %>%
  filter(!ident %in% excluded_briefp_inhibit)
all_dist_briefp_inhibit <- bdd_final_briefp_inhibit %>%
  select(all_of(.[["ident"]])) %>%
  as.dist()

#### filter brief p shift (because of the NA) ----
excluded_briefp_shift <- bdd_final %>% select(ident, ch_briefpshift_y3_cat) %>% filter_all(any_vars(is.na(.)))
excluded_briefp_shift <- excluded_briefp_shift$ident

bdd_final_briefp_shift <- bdd_final %>% 
  select(-all_of(excluded_briefp_shift)) %>%
  filter(!ident %in% excluded_briefp_shift)
all_dist_briefp_shift <- bdd_final_briefp_shift %>%
  select(all_of(.[["ident"]])) %>%
  as.dist()

#### filter brief p emocontrol (because of the NA) ----
excluded_briefp_emocontrol <- bdd_final %>% select(ident, ch_briefpemocontrol_y3_cat) %>% filter_all(any_vars(is.na(.)))
excluded_briefp_emocontrol <- excluded_briefp_emocontrol$ident

bdd_final_briefp_emocontrol <- bdd_final %>% 
  select(-all_of(excluded_briefp_emocontrol)) %>%
  filter(!ident %in% excluded_briefp_emocontrol)
all_dist_briefp_emocontrol <- bdd_final_briefp_emocontrol %>%
  select(all_of(.[["ident"]])) %>%
  as.dist()

#### filter brief p workmemo (because of the NA) ----
excluded_briefp_workmemo <- bdd_final %>% select(ident, ch_briefpworkmemo_y3_cat) %>% filter_all(any_vars(is.na(.)))
excluded_briefp_workmemo <- excluded_briefp_workmemo$ident

bdd_final_briefp_workmemo <- bdd_final %>% 
  select(-all_of(excluded_briefp_workmemo)) %>%
  filter(!ident %in% excluded_briefp_workmemo)
all_dist_briefp_workmemo <- bdd_final_briefp_workmemo %>%
  select(all_of(.[["ident"]])) %>%
  as.dist()

#### filter brief p plan (because of the NA) ----
excluded_briefp_plan <- bdd_final %>% select(ident, ch_briefpplan_y3_cat) %>% filter_all(any_vars(is.na(.)))
excluded_briefp_plan <- excluded_briefp_plan$ident

bdd_final_briefp_plan <- bdd_final %>% 
  select(-all_of(excluded_briefp_plan)) %>%
  filter(!ident %in% excluded_briefp_plan)
all_dist_briefp_plan <- bdd_final_briefp_plan %>%
  select(all_of(.[["ident"]])) %>%
  as.dist()

#### filter IQ (because of the NA) ----
excluded_iq <- bdd_final %>% 
  select(ident, 
         ch_verbal_comprehension_IQ_Y3_cat, 
         ch_visuospatiale_IQ_Y3_cat, 
         ch_work_memory_IQ_Y3_cat, 
         ch_total_IQ_Y3_cat) %>% 
  filter_all(any_vars(is.na(.)))
excluded_iq <- excluded_iq$ident

bdd_final_iq <- bdd_final %>% 
  select(-all_of(excluded_iq)) %>%
  filter(!ident %in% excluded_iq)
all_dist_iq <- bdd_final_iq %>%
  select(all_of(.[["ident"]])) %>%
  as.dist()


## Statistical analysis Adonis2 ----
### CBCL int ----
results_betadiv_univar_bray_curtis_cbcl_int <- 
  adonis2(formula = all_dist_cbcl_int ~ ch_cbclintscore_y2_cat, 
          data = bdd_final_cbcl_int, 
          permutations = 999)
results_betadiv_univar_bray_curtis_cbcl_int <- 
  results_betadiv_univar_bray_curtis_cbcl_int %>%
  rownames_to_column(var = "Variables")%>%
  mutate(Outcome = c(rep("ch_cbclintscore_y2_cat", times = 3)))%>%
  select(Outcome, everything())

results_betadiv_multivar_bray_curtis_cbcl_int <-
  adonis2(formula = reformulate(c("ch_cbclintscore_y2_cat", covariates), response = "all_dist_cbcl_int"),
          data = bdd_final_cbcl_int, 
          permutations = 999)
results_betadiv_multivar_bray_curtis_cbcl_int <- 
  results_betadiv_multivar_bray_curtis_cbcl_int %>%
  rownames_to_column(var = "Variables") %>%
  mutate(Outcome = c(rep("ch_cbclintscore_y2_cat", times = 21))) %>%
  select(Outcome, everything())

### CBCL ext ----
results_betadiv_univar_bray_curtis_cbcl_ext <- 
  adonis2(formula = all_dist_cbcl_ext ~ ch_cbclextscore_y2_cat, 
          data = bdd_final_cbcl_ext, 
          permutations = 999)
results_betadiv_univar_bray_curtis_cbcl_ext <- 
  results_betadiv_univar_bray_curtis_cbcl_ext %>%
  rownames_to_column(var = "Variables")%>%
  mutate(Outcome = c(rep("ch_cbclextscore_y2_cat", times = 3)))%>%
  select(Outcome, everything())

results_betadiv_multivar_bray_curtis_cbcl_ext <-
  adonis2(formula = reformulate(c("ch_cbclextscore_y2_cat", covariates), response = "all_dist_cbcl_ext"),
          data = bdd_final_cbcl_ext, 
          permutations = 999)
results_betadiv_multivar_bray_curtis_cbcl_ext <- 
  results_betadiv_multivar_bray_curtis_cbcl_ext %>%
  rownames_to_column(var = "Variables") %>%
  mutate(Outcome = c(rep("ch_cbclextscore_y2_cat", times = 21))) %>%
  select(Outcome, everything())

### SRS total ----
results_betadiv_univar_bray_curtis_srs <- 
  adonis2(formula = all_dist_srs ~ ch_SRStotal_y3_cat, 
          data = bdd_final_srs, 
          permutations = 999)
results_betadiv_univar_bray_curtis_srs <- 
  results_betadiv_univar_bray_curtis_srs %>%
  rownames_to_column(var = "Variables")%>%
  mutate(Outcome = c(rep("ch_SRStotal_y3_cat", times = 3)))%>%
  select(Outcome, everything())

results_betadiv_multivar_bray_curtis_srs <-
  adonis2(formula = reformulate(c("ch_SRStotal_y3_cat", covariates), response = "all_dist_srs"),
          data = bdd_final_srs, 
          permutations = 999)
results_betadiv_multivar_bray_curtis_srs <- 
  results_betadiv_multivar_bray_curtis_srs %>%
  rownames_to_column(var = "Variables") %>%
  mutate(Outcome = c(rep("ch_SRStotal_y3_cat", times = 21))) %>%
  select(Outcome, everything())

### brief p inhibit ----
results_betadiv_univar_bray_curtis_briefp_inhibit <- 
  adonis2(formula = all_dist_briefp_inhibit ~ ch_briefpinhibit_y3_cat, 
          data = bdd_final_briefp_inhibit, 
          permutations = 999)
results_betadiv_univar_bray_curtis_briefp_inhibit <- 
  results_betadiv_univar_bray_curtis_briefp_inhibit %>%
  rownames_to_column(var = "Variables")%>%
  mutate(Outcome = c(rep("ch_briefpinhibit_y3_cat", times = 3)))%>%
  select(Outcome, everything())

results_betadiv_multivar_bray_curtis_briefp_inhibit <-
  adonis2(formula = reformulate(c("ch_briefpinhibit_y3_cat", covariates), response = "all_dist_briefp_inhibit"),
          data = bdd_final_briefp_inhibit, 
          permutations = 999)
results_betadiv_multivar_bray_curtis_briefp_inhibit <- 
  results_betadiv_multivar_bray_curtis_briefp_inhibit %>%
  rownames_to_column(var = "Variables") %>%
  mutate(Outcome = c(rep("ch_briefpinhibit_y3_cat", times = 21))) %>%
  select(Outcome, everything())

### brief p shift ----
results_betadiv_univar_bray_curtis_briefp_shift <- 
  adonis2(formula = all_dist_briefp_shift ~ ch_briefpshift_y3_cat, 
          data = bdd_final_briefp_shift, 
          permutations = 999)
results_betadiv_univar_bray_curtis_briefp_shift <- 
  results_betadiv_univar_bray_curtis_briefp_shift %>%
  rownames_to_column(var = "Variables")%>%
  mutate(Outcome = c(rep("ch_briefpshift_y3_cat", times = 3)))%>%
  select(Outcome, everything())

results_betadiv_multivar_bray_curtis_briefp_shift <-
  adonis2(formula = reformulate(c("ch_briefpshift_y3_cat", covariates), response = "all_dist_briefp_shift"),
          data = bdd_final_briefp_shift, 
          permutations = 999)
results_betadiv_multivar_bray_curtis_briefp_shift <- 
  results_betadiv_multivar_bray_curtis_briefp_shift %>%
  rownames_to_column(var = "Variables") %>%
  mutate(Outcome = c(rep("ch_briefpshift_y3_cat", times = 21))) %>%
  select(Outcome, everything())

### brief p emocontrol ----
results_betadiv_univar_bray_curtis_briefp_emocontrol <- 
  adonis2(formula = all_dist_briefp_emocontrol ~ ch_briefpemocontrol_y3_cat, 
          data = bdd_final_briefp_emocontrol, 
          permutations = 999)
results_betadiv_univar_bray_curtis_briefp_emocontrol <- 
  results_betadiv_univar_bray_curtis_briefp_emocontrol %>%
  rownames_to_column(var = "Variables")%>%
  mutate(Outcome = c(rep("ch_briefpemocontrol_y3_cat", times = 3)))%>%
  select(Outcome, everything())

results_betadiv_multivar_bray_curtis_briefp_emocontrol <-
  adonis2(formula = reformulate(c("ch_briefpemocontrol_y3_cat", covariates), response = "all_dist_briefp_emocontrol"),
          data = bdd_final_briefp_emocontrol, 
          permutations = 999)
results_betadiv_multivar_bray_curtis_briefp_emocontrol <- 
  results_betadiv_multivar_bray_curtis_briefp_emocontrol %>%
  rownames_to_column(var = "Variables") %>%
  mutate(Outcome = c(rep("ch_briefpemocontrol_y3_cat", times = 21))) %>%
  select(Outcome, everything())

### brief p workmemo ----
results_betadiv_univar_bray_curtis_briefp_workmemo <- 
  adonis2(formula = all_dist_briefp_workmemo ~ ch_briefpworkmemo_y3_cat, 
          data = bdd_final_briefp_workmemo, 
          permutations = 999)
results_betadiv_univar_bray_curtis_briefp_workmemo <- 
  results_betadiv_univar_bray_curtis_briefp_workmemo %>%
  rownames_to_column(var = "Variables")%>%
  mutate(Outcome = c(rep("ch_briefpworkmemo_y3_cat", times = 3)))%>%
  select(Outcome, everything())

results_betadiv_multivar_bray_curtis_briefp_workmemo <-
  adonis2(formula = reformulate(c("ch_briefpworkmemo_y3_cat", covariates), response = "all_dist_briefp_workmemo"),
          data = bdd_final_briefp_workmemo, 
          permutations = 999)
results_betadiv_multivar_bray_curtis_briefp_workmemo <- 
  results_betadiv_multivar_bray_curtis_briefp_workmemo %>%
  rownames_to_column(var = "Variables") %>%
  mutate(Outcome = c(rep("ch_briefpworkmemo_y3_cat", times = 21))) %>%
  select(Outcome, everything())

### brief p plan ----
results_betadiv_univar_bray_curtis_briefp_plan <- 
  adonis2(formula = all_dist_briefp_plan ~ ch_briefpplan_y3_cat, 
          data = bdd_final_briefp_plan, 
          permutations = 999)
results_betadiv_univar_bray_curtis_briefp_plan <- 
  results_betadiv_univar_bray_curtis_briefp_plan %>%
  rownames_to_column(var = "Variables")%>%
  mutate(Outcome = c(rep("ch_briefpplan_y3_cat", times = 3)))%>%
  select(Outcome, everything())

results_betadiv_multivar_bray_curtis_briefp_plan <-
  adonis2(formula = reformulate(c("ch_briefpplan_y3_cat", covariates), response = "all_dist_briefp_plan"),
          data = bdd_final_briefp_plan, 
          permutations = 999)
results_betadiv_multivar_bray_curtis_briefp_plan <- 
  results_betadiv_multivar_bray_curtis_briefp_plan %>%
  rownames_to_column(var = "Variables") %>%
  mutate(Outcome = c(rep("ch_briefpplan_y3_cat", times = 21))) %>%
  select(Outcome, everything())

### IQ ----
iq_vars <- 
  bdd_final_iq%>% 
  select(ch_verbal_comprehension_IQ_Y3_cat, 
         ch_visuospatiale_IQ_Y3_cat, 
         ch_work_memory_IQ_Y3_cat, 
         ch_total_IQ_Y3_cat) %>%
  colnames()

results_betadiv_univar_bray_curtis_iq <- 
  lapply(iq_vars, function(x) {
    formula <- reformulate(x, response = "all_dist_iq")
    adonis2(formula, data = bdd_final_iq, permutations = 999)
  })
results_betadiv_univar_bray_curtis_iq <- 
  do.call(rbind, results_betadiv_univar_bray_curtis_iq) %>%
  rownames_to_column(var = "Variables") %>%
  mutate(
    Outcome = c(rep("ch_verbal_comprehension_IQ_Y3_cat", times = 3), 
                   rep("ch_visuospatiale_IQ_Y3_cat", times = 3), 
                   rep("ch_work_memory_IQ_Y3_cat", times = 3),
                   rep("ch_total_IQ_Y3_cat", times = 3))) %>%
  select(Outcome, everything())

results_betadiv_multivar_bray_curtis_iq <-
  lapply(iq_vars, function(x) {
    formula <- reformulate(c(x, covariates), response = "all_dist_iq")
    adonis2(formula, data = bdd_final_iq, permutations = 999)
  })
results_betadiv_multivar_bray_curtis_iq <-
  do.call(rbind, results_betadiv_multivar_bray_curtis_iq) %>%
  rownames_to_column(var = "Variables") %>%
  mutate(
    Outcome = c(rep("ch_verbal_comprehension_IQ_Y3_cat", times = 21), 
                 rep("ch_visuospatiale_IQ_Y3_cat", times = 21), 
                 rep("ch_work_memory_IQ_Y3_cat", times = 21),
                 rep("ch_total_IQ_Y3_cat", times = 21))) %>%
  select(Outcome, everything())

results_betadiv_univar_bray_curtis_iq %>%                         # Visualisation des résultats significatifs univarié
  filter(`Pr(>F)` < 0.05) %>%
  filter(`Variables` %in% c(iq_vars)) %>%
  View()

results_betadiv_multivar_bray_curtis_iq %>%                     # Visualisation des résultats significatifs multivarié
  filter(`Pr(>F)` < 0.05) %>%
  filter(`Variables` %in% c(iq_vars)) %>%
  View()



### Assemblage ----
results_betadiv_univ <-
  list(
    results_betadiv_univar_bray_curtis_cbcl_int,
    results_betadiv_univar_bray_curtis_cbcl_ext,
    results_betadiv_univar_bray_curtis_srs, 
    results_betadiv_univar_bray_curtis_briefp_inhibit, 
    results_betadiv_univar_bray_curtis_briefp_shift, 
    results_betadiv_univar_bray_curtis_briefp_emocontrol, 
    results_betadiv_univar_bray_curtis_briefp_workmemo, 
    results_betadiv_univar_bray_curtis_briefp_plan, 
    results_betadiv_univar_bray_curtis_iq)
results_betadiv_univ <- do.call(rbind, results_betadiv_univ, quote = FALSE)
results_betadiv_univ <- results_betadiv_univ %>%
  mutate(Variables = if_else(Variables %in% outcomes, Variables, str_remove(Variables, "\\d+$")))

results_betadiv_multi <-
  list(
    results_betadiv_multivar_bray_curtis_cbcl_int,
    results_betadiv_multivar_bray_curtis_cbcl_ext,
    results_betadiv_multivar_bray_curtis_srs, 
    results_betadiv_multivar_bray_curtis_briefp_inhibit, 
    results_betadiv_multivar_bray_curtis_briefp_shift, 
    results_betadiv_multivar_bray_curtis_briefp_emocontrol, 
    results_betadiv_multivar_bray_curtis_briefp_workmemo, 
    results_betadiv_multivar_bray_curtis_briefp_plan, 
    results_betadiv_multivar_bray_curtis_iq)
results_betadiv_multi <- do.call(rbind, results_betadiv_multi, quote = FALSE)
results_betadiv_multi <- results_betadiv_multi %>%
  mutate(Variables = if_else(Variables %in% outcomes, Variables, str_remove(Variables, "\\d+$"))) %>%
  filter(Variables %in% outcomes | str_detect(Variables, "Residual") | str_detect(Variables, "Total"))


results_betadiv_complet <-
  list(
    univar = list(results_betadiv_univar_bray_curtis_cbcl_int = results_betadiv_univar_bray_curtis_cbcl_int,
                  results_betadiv_univar_bray_curtis_cbcl_ext = results_betadiv_univar_bray_curtis_cbcl_ext,
                  results_betadiv_univar_bray_curtis_srs = results_betadiv_univar_bray_curtis_srs, 
                  results_betadiv_univar_bray_curtis_briefp_inhibit = results_betadiv_univar_bray_curtis_briefp_inhibit, 
                  results_betadiv_univar_bray_curtis_briefp_shift = results_betadiv_univar_bray_curtis_briefp_shift, 
                  results_betadiv_univar_bray_curtis_briefp_emocontrol = results_betadiv_univar_bray_curtis_briefp_emocontrol, 
                  results_betadiv_univar_bray_curtis_briefp_workmemo =results_betadiv_univar_bray_curtis_briefp_workmemo, 
                  results_betadiv_univar_bray_curtis_briefp_plan = results_betadiv_univar_bray_curtis_briefp_plan, 
                  results_betadiv_univar_bray_curtis_iq = results_betadiv_univar_bray_curtis_iq),
    multivar = list(results_betadiv_multivar_bray_curtis_cbcl_int = results_betadiv_multivar_bray_curtis_cbcl_int,
                    results_betadiv_multivar_bray_curtis_cbcl_ext = results_betadiv_multivar_bray_curtis_cbcl_ext,
                    results_betadiv_multivar_bray_curtis_srs = results_betadiv_multivar_bray_curtis_srs, 
                    results_betadiv_multivar_bray_curtis_briefp_inhibit = results_betadiv_multivar_bray_curtis_briefp_inhibit, 
                    results_betadiv_multivar_bray_curtis_briefp_shift = results_betadiv_multivar_bray_curtis_briefp_shift, 
                    results_betadiv_multivar_bray_curtis_briefp_emocontrol = results_betadiv_multivar_bray_curtis_briefp_emocontrol, 
                    results_betadiv_multivar_bray_curtis_briefp_workmemo = results_betadiv_multivar_bray_curtis_briefp_workmemo, 
                    results_betadiv_multivar_bray_curtis_briefp_plan = results_betadiv_multivar_bray_curtis_briefp_plan, 
                    results_betadiv_multivar_bray_curtis_iq = results_betadiv_multivar_bray_curtis_iq))

write_xlsx(list(results_betadiv_multi = results_betadiv_multi, 
                results_betadiv_univ = results_betadiv_univ), 
           path = "4_output/results_betadiv.xlsx")

rm(excluded_cbcl_int, 
   excluded_cbcl_ext, 
   excluded_srs, 
   excluded_briefp_inhibit, 
   excluded_briefp_shift, 
   excluded_briefp_emocontrol, 
   excluded_briefp_workmemo, 
   excluded_briefp_plan, 
   excluded_iq, 
   
   bdd_final_cbcl_int, 
   bdd_final_cbcl_ext, 
   bdd_final_srs, 
   bdd_final_briefp_inhibit, 
   bdd_final_briefp_shift, 
   bdd_final_briefp_emocontrol, 
   bdd_final_briefp_workmemo, 
   bdd_final_briefp_plan, 
   bdd_final_iq, 
   
   all_dist_cbcl_int, 
   all_dist_cbcl_ext, 
   all_dist_srs, 
   all_dist_briefp_inhibit, 
   all_dist_briefp_shift, 
   all_dist_briefp_emocontrol, 
   all_dist_briefp_workmemo, 
   all_dist_briefp_plan, 
   all_dist_iq, 
   
   iq_vars,
   nmds,
   
   results_betadiv_univar_bray_curtis_cbcl_int,
   results_betadiv_univar_bray_curtis_cbcl_ext,
   results_betadiv_univar_bray_curtis_srs, 
   results_betadiv_univar_bray_curtis_briefp_inhibit, 
   results_betadiv_univar_bray_curtis_briefp_shift, 
   results_betadiv_univar_bray_curtis_briefp_emocontrol, 
   results_betadiv_univar_bray_curtis_briefp_workmemo, 
   results_betadiv_univar_bray_curtis_briefp_plan, 
   results_betadiv_univar_bray_curtis_iq, 
   
   results_betadiv_multivar_bray_curtis_cbcl_int,
   results_betadiv_multivar_bray_curtis_cbcl_ext,
   results_betadiv_multivar_bray_curtis_srs, 
   results_betadiv_multivar_bray_curtis_briefp_inhibit, 
   results_betadiv_multivar_bray_curtis_briefp_shift, 
   results_betadiv_multivar_bray_curtis_briefp_emocontrol, 
   results_betadiv_multivar_bray_curtis_briefp_workmemo, 
   results_betadiv_multivar_bray_curtis_briefp_plan, 
   results_betadiv_multivar_bray_curtis_iq)


results_betadiv_univ %>%
  filter(`Pr(>F)` < 0.05) %>%
  filter(Variables %in% outcomes) %>%
  View()

results_betadiv_multi %>%
  filter(`Pr(>F)` < 0.05) %>%
  filter(Variables %in% outcomes) %>%
  View()

table_3 <- results_betadiv_multi %>%
  filter(!Variables %in% c("Residual", "Total")) %>%
  select(-Variables, -Df, -R2) %>%
  mutate(
    Outcome = fct_recode(Outcome, 
                         "Internalizing CBCL score at 2 years" = "ch_cbclintscore_y2_cat", 
                         "Externalizing CBCL score at 2 years" = "ch_cbclextscore_y2_cat", 
                         "Total SRS score at 3 years" = "ch_SRStotal_y3_cat", 
                         "Inhibition BRIEF-P score at 3 years" = "ch_briefpinhibit_y3_cat", 
                         "Shift BRIEF-P score at 3 years" = "ch_briefpshift_y3_cat", 
                         "Emotional control BRIEF-P score at 3 years" = "ch_briefpemocontrol_y3_cat",
                         "Working memory BRIEF-P score at 3 years" = "ch_briefpworkmemo_y3_cat",
                         "Plan and organization BRIEF-P score at 3 years" = "ch_briefpplan_y3_cat",
                         
                         "Verbal comprehension IQ score at 3 years" = "ch_verbal_comprehension_IQ_Y3_cat", 
                         "Visuospatial IQ score at 3 years" = "ch_visuospatiale_IQ_Y3_cat",
                         "Work memory IQ score at 3 years" = "ch_work_memory_IQ_Y3_cat",
                         "Total IQ score at 3 years" = "ch_total_IQ_Y3_cat"))
write_xlsx(table_3, 
           path = "4_output/table_3.xlsx")

## Boxplots ----
### Trim2. ----
# on réduit la matrice de dissimilarité pour ne pas avoir les données en double
metric_bray_curtis_red <- metric_bray_curtis %>% column_to_rownames(var = "ident")
metric_bray_curtis_red[lower.tri(metric_bray_curtis_red)]<- NA
metric_bray_curtis_red <- metric_bray_curtis_red %>% rownames_to_column(var = "ident")

# on merge les données d'expo t2 aux données de dissimilarité
bdd_long_t2 <- bdd_final_t2 %>% select(ident, all_of(explanatory_vars_t2))
bdd_long_t2 <- full_join(bdd_long_t2, metric_bray_curtis_red, by = "ident")

# on fait passer la base de donnnées en long
bdd_long_t2[,explanatory_vars_t2] <- lapply(bdd_long_t2[,explanatory_vars_t2], as.character)
bdd_long_t2 <- bdd_long_t2 %>%
  pivot_longer(cols = c(-"ident", -all_of(explanatory_vars_t2)), 
               names_to = "ident_bis", 
               values_to = "Bray_curtis_dissimilarity") %>%
  select(ident, ident_bis, Bray_curtis_dissimilarity, everything()) %>%
  filter(ident != ident_bis) %>%
  filter(!is.na(Bray_curtis_dissimilarity)) 

# on créé une base de données bis pour savoir à quelle groupes d'expos sont comparés les paires 
bdd_long_t2_bis <- bdd_final_t2 %>%
  select(ident, 
         all_of(explanatory_vars_t2)) %>%
  rename(ident_bis = ident, 
         mo_DEHP_ms_i_cor_t2_ter_bis = mo_DEHP_ms_i_cor_t2_ter, 
         mo_MnBP_i_cor_t2_ter_bis = mo_MnBP_i_cor_t2_ter, 
         mo_DiNP_ms_i_cor_t2_ter_bis = mo_DiNP_ms_i_cor_t2_ter, 
         mo_MiBP_i_cor_t2_ter_bis = mo_MiBP_i_cor_t2_ter, 
         mo_MBzP_i_cor_t2_ter_bis = mo_MBzP_i_cor_t2_ter, 
         mo_MEP_i_cor_t2_ter_bis = mo_MEP_i_cor_t2_ter, 
         mo_ohMPHP_i_cor_t2_ter_bis = mo_ohMPHP_i_cor_t2_ter, 
         mo_DINCH_ms_i_cor_t2_ter_bis = mo_DINCH_ms_i_cor_t2_ter)

# on merge les données puis on supprime la base de données bis 
bdd_long_t2 <- full_join(bdd_long_t2, bdd_long_t2_bis, by = "ident_bis") 
bdd_long_t2 <- bdd_long_t2 %>%
  filter(ident != 15804) %>%
  filter(ident_bis != 15804)
rm(bdd_long_t2_bis)

# on créé une variable qui indique de quel intergroupe il s'agit
bdd_long_t2 <- bdd_long_t2 %>%
  mutate(
    DEHP = case_when(mo_DEHP_ms_i_cor_t2_ter == "1st tertile" & mo_DEHP_ms_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile", 
                     mo_DEHP_ms_i_cor_t2_ter == "2nd tertile" & mo_DEHP_ms_i_cor_t2_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     mo_DEHP_ms_i_cor_t2_ter == "3rd tertile" & mo_DEHP_ms_i_cor_t2_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     mo_DEHP_ms_i_cor_t2_ter == "1st tertile" & mo_DEHP_ms_i_cor_t2_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     mo_DEHP_ms_i_cor_t2_ter == "2nd tertile" & mo_DEHP_ms_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     mo_DEHP_ms_i_cor_t2_ter == "1st tertile" & mo_DEHP_ms_i_cor_t2_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     mo_DEHP_ms_i_cor_t2_ter == "3rd tertile" & mo_DEHP_ms_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     mo_DEHP_ms_i_cor_t2_ter == "2nd tertile" & mo_DEHP_ms_i_cor_t2_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     mo_DEHP_ms_i_cor_t2_ter == "3rd tertile" & mo_DEHP_ms_i_cor_t2_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    MnBP = case_when(mo_MnBP_i_cor_t2_ter == "1st tertile" & mo_MnBP_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile", 
                     mo_MnBP_i_cor_t2_ter == "2nd tertile" & mo_MnBP_i_cor_t2_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     mo_MnBP_i_cor_t2_ter == "3rd tertile" & mo_MnBP_i_cor_t2_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     mo_MnBP_i_cor_t2_ter == "1st tertile" & mo_MnBP_i_cor_t2_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     mo_MnBP_i_cor_t2_ter == "2nd tertile" & mo_MnBP_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     mo_MnBP_i_cor_t2_ter == "1st tertile" & mo_MnBP_i_cor_t2_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     mo_MnBP_i_cor_t2_ter == "3rd tertile" & mo_MnBP_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     mo_MnBP_i_cor_t2_ter == "2nd tertile" & mo_MnBP_i_cor_t2_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     mo_MnBP_i_cor_t2_ter == "3rd tertile" & mo_MnBP_i_cor_t2_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    DiNP = case_when(mo_DiNP_ms_i_cor_t2_ter == "1st tertile" & mo_DiNP_ms_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile", 
                     mo_DiNP_ms_i_cor_t2_ter == "2nd tertile" & mo_DiNP_ms_i_cor_t2_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     mo_DiNP_ms_i_cor_t2_ter == "3rd tertile" & mo_DiNP_ms_i_cor_t2_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     mo_DiNP_ms_i_cor_t2_ter == "1st tertile" & mo_DiNP_ms_i_cor_t2_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     mo_DiNP_ms_i_cor_t2_ter == "2nd tertile" & mo_DiNP_ms_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     mo_DiNP_ms_i_cor_t2_ter == "1st tertile" & mo_DiNP_ms_i_cor_t2_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     mo_DiNP_ms_i_cor_t2_ter == "3rd tertile" & mo_DiNP_ms_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     mo_DiNP_ms_i_cor_t2_ter == "2nd tertile" & mo_DiNP_ms_i_cor_t2_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     mo_DiNP_ms_i_cor_t2_ter == "3rd tertile" & mo_DiNP_ms_i_cor_t2_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    MiBP = case_when(mo_MiBP_i_cor_t2_ter == "1st tertile" & mo_MiBP_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile", 
                     mo_MiBP_i_cor_t2_ter == "2nd tertile" & mo_MiBP_i_cor_t2_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     mo_MiBP_i_cor_t2_ter == "3rd tertile" & mo_MiBP_i_cor_t2_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     mo_MiBP_i_cor_t2_ter == "1st tertile" & mo_MiBP_i_cor_t2_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     mo_MiBP_i_cor_t2_ter == "2nd tertile" & mo_MiBP_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     mo_MiBP_i_cor_t2_ter == "1st tertile" & mo_MiBP_i_cor_t2_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     mo_MiBP_i_cor_t2_ter == "3rd tertile" & mo_MiBP_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     mo_MiBP_i_cor_t2_ter == "2nd tertile" & mo_MiBP_i_cor_t2_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     mo_MiBP_i_cor_t2_ter == "3rd tertile" & mo_MiBP_i_cor_t2_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    MBzP = case_when(mo_MBzP_i_cor_t2_ter == "1st tertile" & mo_MBzP_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile", 
                     mo_MBzP_i_cor_t2_ter == "2nd tertile" & mo_MBzP_i_cor_t2_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     mo_MBzP_i_cor_t2_ter == "3rd tertile" & mo_MBzP_i_cor_t2_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     mo_MBzP_i_cor_t2_ter == "1st tertile" & mo_MBzP_i_cor_t2_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     mo_MBzP_i_cor_t2_ter == "2nd tertile" & mo_MBzP_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     mo_MBzP_i_cor_t2_ter == "1st tertile" & mo_MBzP_i_cor_t2_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     mo_MBzP_i_cor_t2_ter == "3rd tertile" & mo_MBzP_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     mo_MBzP_i_cor_t2_ter == "2nd tertile" & mo_MBzP_i_cor_t2_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     mo_MBzP_i_cor_t2_ter == "3rd tertile" & mo_MBzP_i_cor_t2_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    MEP = case_when(mo_MEP_i_cor_t2_ter == "1st tertile" & mo_MEP_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile", 
                    mo_MEP_i_cor_t2_ter == "2nd tertile" & mo_MEP_i_cor_t2_ter_bis == "2nd tertile" ~ "2nd tertile", 
                    mo_MEP_i_cor_t2_ter == "3rd tertile" & mo_MEP_i_cor_t2_ter_bis == "3rd tertile" ~ "3rd tertile", 
                    
                    mo_MEP_i_cor_t2_ter == "1st tertile" & mo_MEP_i_cor_t2_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                    mo_MEP_i_cor_t2_ter == "2nd tertile" & mo_MEP_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                    
                    mo_MEP_i_cor_t2_ter == "1st tertile" & mo_MEP_i_cor_t2_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                    mo_MEP_i_cor_t2_ter == "3rd tertile" & mo_MEP_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                    
                    mo_MEP_i_cor_t2_ter == "2nd tertile" & mo_MEP_i_cor_t2_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                    mo_MEP_i_cor_t2_ter == "3rd tertile" & mo_MEP_i_cor_t2_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    ohMPHP = case_when(mo_ohMPHP_i_cor_t2_ter == "1st tertile" & mo_ohMPHP_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile", 
                       mo_ohMPHP_i_cor_t2_ter == "2nd tertile" & mo_ohMPHP_i_cor_t2_ter_bis == "2nd tertile" ~ "2nd tertile", 
                       mo_ohMPHP_i_cor_t2_ter == "3rd tertile" & mo_ohMPHP_i_cor_t2_ter_bis == "3rd tertile" ~ "3rd tertile", 
                       
                       mo_ohMPHP_i_cor_t2_ter == "1st tertile" & mo_ohMPHP_i_cor_t2_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                       mo_ohMPHP_i_cor_t2_ter == "2nd tertile" & mo_ohMPHP_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                       
                       mo_ohMPHP_i_cor_t2_ter == "1st tertile" & mo_ohMPHP_i_cor_t2_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                       mo_ohMPHP_i_cor_t2_ter == "3rd tertile" & mo_ohMPHP_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                       
                       mo_ohMPHP_i_cor_t2_ter == "2nd tertile" & mo_ohMPHP_i_cor_t2_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                       mo_ohMPHP_i_cor_t2_ter == "3rd tertile" & mo_ohMPHP_i_cor_t2_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    DINCH = case_when(mo_DINCH_ms_i_cor_t2_ter == "1st tertile" & mo_DINCH_ms_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile", 
                     mo_DINCH_ms_i_cor_t2_ter == "2nd tertile" & mo_DINCH_ms_i_cor_t2_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     mo_DINCH_ms_i_cor_t2_ter == "3rd tertile" & mo_DINCH_ms_i_cor_t2_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     mo_DINCH_ms_i_cor_t2_ter == "1st tertile" & mo_DINCH_ms_i_cor_t2_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     mo_DINCH_ms_i_cor_t2_ter == "2nd tertile" & mo_DINCH_ms_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     mo_DINCH_ms_i_cor_t2_ter == "1st tertile" & mo_DINCH_ms_i_cor_t2_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     mo_DINCH_ms_i_cor_t2_ter == "3rd tertile" & mo_DINCH_ms_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     mo_DINCH_ms_i_cor_t2_ter == "2nd tertile" & mo_DINCH_ms_i_cor_t2_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     mo_DINCH_ms_i_cor_t2_ter == "3rd tertile" & mo_DINCH_ms_i_cor_t2_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"))%>%
  select(ident, ident_bis, Bray_curtis_dissimilarity, 
         DEHP, MnBP, DiNP, MiBP, MBzP, MEP, ohMPHP, DINCH)



# on fait passer les données en long 
bdd_long_t2 <- bdd_long_t2 %>%
  pivot_longer(cols = -c("ident", "ident_bis", "Bray_curtis_dissimilarity"), 
               names_to = "Pollutant", 
               values_to = "Groups")


# duplication de 1st tertile, 2nd tertile, 3rd tertile
test <- 
  bdd_long_t2 %>% 
  filter(Groups %in% c("1st tertile", "2nd tertile", "3rd tertile")) %>%  
  mutate(Groups = fct_recode(Groups, 
                             " 1st tertile" = "1st tertile", 
                             " 2nd tertile" = "2nd tertile", 
                             " 3rd tertile" = "3rd tertile"))
bdd_long_t2 <- bind_rows(bdd_long_t2, test)
rm(test)

bdd_long_t2 <- bdd_long_t2 %>%
  mutate(
    Pollutant_rec = fct_recode(Pollutant, 
                               "ΣDEHP trim.2" = "DEHP",
                                "ΣDINCH trim.2" = "DINCH",
                                "ΣDiNP trim.2" = "DiNP",
                                "MBzP trim.2" = "MBzP",
                                "MEP trim.2" = "MEP",
                                "MiBP trim.2" = "MiBP",
                                "MnBP trim.2" = "MnBP",
                                "ohMPHP trim.2" = "ohMPHP"),
    Pollutant = fct_relevel(Pollutant,
                            "DEHP", "MnBP", "DiNP", "MiBP", "MBzP", "MEP", "ohMPHP", "DINCH"),
    Pollutant_rec = fct_relevel(Pollutant_rec,
                            "ΣDEHP trim.2", "MnBP trim.2", "ΣDiNP trim.2", "MiBP trim.2", 
                            "MBzP trim.2", "MEP trim.2", "ohMPHP trim.2", "ΣDINCH trim.2"),
    Groups = fct_relevel(Groups,
                         " 3rd tertile",                         # intra groupe forte exposition 
                         "2nd tertile - 3rd tertile",            # inter groupe moyenne et forte exposition
                         " 2nd tertile",                         # intra groupe moyenne exposition 
                         
                         "3rd tertile",                          # intra groupe forte exposition 
                         "1st tertile - 3rd tertile",            # inter groupe faible et forte exposition 
                         " 1st tertile", 
                         
                         "2nd tertile",                          # intra groupe moyenne exposition 
                         "1st tertile - 2nd tertile",            # inter groupe faible et moyenne exposition 
                         "1st tertile"),                         # intra groupe faible exposition
    Groups_rec = fct_recode(Groups, 
                            "Medium vs High" = " 3rd tertile",
                            "Medium vs High" = "2nd tertile - 3rd tertile",
                            "Medium vs High" = " 2nd tertile",
                            "Low vs High" = "3rd tertile",
                            "Low vs High" = "1st tertile - 3rd tertile",
                            "Low vs High" = " 1st tertile",
                            "Low vs Medium" = "2nd tertile",
                            "Low vs Medium" = "1st tertile - 2nd tertile",
                            "Low vs Medium" = "1st tertile"))



#### boxplot version 1 ----
boxplot_t2_phthalates <- bdd_long_t2 %>%                    # intra groupe faible exposition
  ggplot() +
  aes(
    x = Bray_curtis_dissimilarity,
    y = Groups,
    fill = Groups
  ) +
  geom_boxplot() +  
  scale_fill_manual(
    values = c(`1st tertile` = "#FFE0DE",   # les catégories faibles expositions (intragroupes)
               ` 1st tertile` = "#FFE0DE", 
               
               `1st tertile - 2nd tertile` = "gray70",  # les variances intergroupes
               `1st tertile - 3rd tertile` = "gray70",
               `2nd tertile - 3rd tertile` = "gray70",
               
               
               `2nd tertile` = "#FF8F87",   # les catégorires moyennes expositions (intragroupes)
               ` 2nd tertile` = "#FF8F87",  
               
               `3rd tertile` = "#FF4034",   
               ` 3rd tertile` = "#FF4034")
  ) +
  labs(x = "Bray Curtis dissimilarity", 
       y = "") +
  theme_lucid() +
  theme(legend.position = "none", 
        axis.title = element_text(face = "bold")) +
  facet_wrap(~Pollutant_rec, scales = "free", ncol = 2)

### Trim3. ----
# on réduit la matrice de dissimilarité pour ne pas avoir les données en double
metric_bray_curtis_red <- metric_bray_curtis %>% column_to_rownames(var = "ident")
metric_bray_curtis_red[lower.tri(metric_bray_curtis_red)]<- NA
metric_bray_curtis_red <- metric_bray_curtis_red %>% rownames_to_column(var = "ident")

# on merge les données d'expo t3 aux données de dissimilarité
bdd_long_t3 <- bdd_final_t3 %>% select(ident, all_of(explanatory_vars_t3))
bdd_long_t3 <- full_join(bdd_long_t3, metric_bray_curtis_red, by = "ident")

# on fait passer la base de donnnées en long
bdd_long_t3[,explanatory_vars_t3] <- lapply(bdd_long_t3[,explanatory_vars_t3], as.character)
bdd_long_t3 <- bdd_long_t3 %>%
  pivot_longer(cols = c(-"ident", -all_of(explanatory_vars_t3)), 
               names_to = "ident_bis", 
               values_to = "Bray_curtis_dissimilarity") %>%
  select(ident, ident_bis, Bray_curtis_dissimilarity, everything()) %>%
  filter(ident != ident_bis) %>%
  filter(!is.na(Bray_curtis_dissimilarity)) 

# on créé une base de données bis pour savoir à quelle groupes d'expos sont comparés les paires 
bdd_long_t3_bis <- bdd_final_t3 %>%
  select(ident, 
         all_of(explanatory_vars_t3)) %>%
  rename(ident_bis = ident, 
         mo_DEHP_ms_i_cor_t3_ter_bis = mo_DEHP_ms_i_cor_t3_ter, 
         mo_MnBP_i_cor_t3_ter_bis = mo_MnBP_i_cor_t3_ter, 
         mo_DiNP_ms_i_cor_t3_ter_bis = mo_DiNP_ms_i_cor_t3_ter, 
         mo_MiBP_i_cor_t3_ter_bis = mo_MiBP_i_cor_t3_ter, 
         mo_MBzP_i_cor_t3_ter_bis = mo_MBzP_i_cor_t3_ter, 
         mo_MEP_i_cor_t3_ter_bis = mo_MEP_i_cor_t3_ter, 
         mo_ohMPHP_i_cor_t3_ter_bis = mo_ohMPHP_i_cor_t3_ter, 
         mo_DINCH_ms_i_cor_t3_ter_bis = mo_DINCH_ms_i_cor_t3_ter)

# on merge les données puis on supprime la base de données bis 
bdd_long_t3 <- full_join(bdd_long_t3, bdd_long_t3_bis, by = "ident_bis") 
bdd_long_t3 <- bdd_long_t3 %>%
  filter(!ident %in% c(17827, 15929, 26891, 25668, 23330, 28199)) %>%
  filter(!ident_bis %in% c(17827, 15929, 26891, 25668, 23330, 28199)) %>%
  filter(!is.na(ident))
rm(bdd_long_t3_bis)

# on créé une variable qui indique de quel intergroupe il s'agit
bdd_long_t3 <- bdd_long_t3 %>%
  mutate(
    DEHP = case_when(mo_DEHP_ms_i_cor_t3_ter == "1st tertile" & mo_DEHP_ms_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile", 
                     mo_DEHP_ms_i_cor_t3_ter == "2nd tertile" & mo_DEHP_ms_i_cor_t3_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     mo_DEHP_ms_i_cor_t3_ter == "3rd tertile" & mo_DEHP_ms_i_cor_t3_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     mo_DEHP_ms_i_cor_t3_ter == "1st tertile" & mo_DEHP_ms_i_cor_t3_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     mo_DEHP_ms_i_cor_t3_ter == "2nd tertile" & mo_DEHP_ms_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     mo_DEHP_ms_i_cor_t3_ter == "1st tertile" & mo_DEHP_ms_i_cor_t3_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     mo_DEHP_ms_i_cor_t3_ter == "3rd tertile" & mo_DEHP_ms_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     mo_DEHP_ms_i_cor_t3_ter == "2nd tertile" & mo_DEHP_ms_i_cor_t3_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     mo_DEHP_ms_i_cor_t3_ter == "3rd tertile" & mo_DEHP_ms_i_cor_t3_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    MnBP = case_when(mo_MnBP_i_cor_t3_ter == "1st tertile" & mo_MnBP_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile", 
                     mo_MnBP_i_cor_t3_ter == "2nd tertile" & mo_MnBP_i_cor_t3_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     mo_MnBP_i_cor_t3_ter == "3rd tertile" & mo_MnBP_i_cor_t3_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     mo_MnBP_i_cor_t3_ter == "1st tertile" & mo_MnBP_i_cor_t3_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     mo_MnBP_i_cor_t3_ter == "2nd tertile" & mo_MnBP_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     mo_MnBP_i_cor_t3_ter == "1st tertile" & mo_MnBP_i_cor_t3_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     mo_MnBP_i_cor_t3_ter == "3rd tertile" & mo_MnBP_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     mo_MnBP_i_cor_t3_ter == "2nd tertile" & mo_MnBP_i_cor_t3_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     mo_MnBP_i_cor_t3_ter == "3rd tertile" & mo_MnBP_i_cor_t3_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    DiNP = case_when(mo_DiNP_ms_i_cor_t3_ter == "1st tertile" & mo_DiNP_ms_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile", 
                     mo_DiNP_ms_i_cor_t3_ter == "2nd tertile" & mo_DiNP_ms_i_cor_t3_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     mo_DiNP_ms_i_cor_t3_ter == "3rd tertile" & mo_DiNP_ms_i_cor_t3_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     mo_DiNP_ms_i_cor_t3_ter == "1st tertile" & mo_DiNP_ms_i_cor_t3_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     mo_DiNP_ms_i_cor_t3_ter == "2nd tertile" & mo_DiNP_ms_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     mo_DiNP_ms_i_cor_t3_ter == "1st tertile" & mo_DiNP_ms_i_cor_t3_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     mo_DiNP_ms_i_cor_t3_ter == "3rd tertile" & mo_DiNP_ms_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     mo_DiNP_ms_i_cor_t3_ter == "2nd tertile" & mo_DiNP_ms_i_cor_t3_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     mo_DiNP_ms_i_cor_t3_ter == "3rd tertile" & mo_DiNP_ms_i_cor_t3_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    MiBP = case_when(mo_MiBP_i_cor_t3_ter == "1st tertile" & mo_MiBP_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile", 
                     mo_MiBP_i_cor_t3_ter == "2nd tertile" & mo_MiBP_i_cor_t3_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     mo_MiBP_i_cor_t3_ter == "3rd tertile" & mo_MiBP_i_cor_t3_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     mo_MiBP_i_cor_t3_ter == "1st tertile" & mo_MiBP_i_cor_t3_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     mo_MiBP_i_cor_t3_ter == "2nd tertile" & mo_MiBP_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     mo_MiBP_i_cor_t3_ter == "1st tertile" & mo_MiBP_i_cor_t3_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     mo_MiBP_i_cor_t3_ter == "3rd tertile" & mo_MiBP_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     mo_MiBP_i_cor_t3_ter == "2nd tertile" & mo_MiBP_i_cor_t3_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     mo_MiBP_i_cor_t3_ter == "3rd tertile" & mo_MiBP_i_cor_t3_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    MBzP = case_when(mo_MBzP_i_cor_t3_ter == "1st tertile" & mo_MBzP_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile", 
                     mo_MBzP_i_cor_t3_ter == "2nd tertile" & mo_MBzP_i_cor_t3_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     mo_MBzP_i_cor_t3_ter == "3rd tertile" & mo_MBzP_i_cor_t3_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     mo_MBzP_i_cor_t3_ter == "1st tertile" & mo_MBzP_i_cor_t3_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     mo_MBzP_i_cor_t3_ter == "2nd tertile" & mo_MBzP_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     mo_MBzP_i_cor_t3_ter == "1st tertile" & mo_MBzP_i_cor_t3_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     mo_MBzP_i_cor_t3_ter == "3rd tertile" & mo_MBzP_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     mo_MBzP_i_cor_t3_ter == "2nd tertile" & mo_MBzP_i_cor_t3_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     mo_MBzP_i_cor_t3_ter == "3rd tertile" & mo_MBzP_i_cor_t3_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    MEP = case_when(mo_MEP_i_cor_t3_ter == "1st tertile" & mo_MEP_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile", 
                    mo_MEP_i_cor_t3_ter == "2nd tertile" & mo_MEP_i_cor_t3_ter_bis == "2nd tertile" ~ "2nd tertile", 
                    mo_MEP_i_cor_t3_ter == "3rd tertile" & mo_MEP_i_cor_t3_ter_bis == "3rd tertile" ~ "3rd tertile", 
                    
                    mo_MEP_i_cor_t3_ter == "1st tertile" & mo_MEP_i_cor_t3_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                    mo_MEP_i_cor_t3_ter == "2nd tertile" & mo_MEP_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                    
                    mo_MEP_i_cor_t3_ter == "1st tertile" & mo_MEP_i_cor_t3_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                    mo_MEP_i_cor_t3_ter == "3rd tertile" & mo_MEP_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                    
                    mo_MEP_i_cor_t3_ter == "2nd tertile" & mo_MEP_i_cor_t3_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                    mo_MEP_i_cor_t3_ter == "3rd tertile" & mo_MEP_i_cor_t3_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    ohMPHP = case_when(mo_ohMPHP_i_cor_t3_ter == "1st tertile" & mo_ohMPHP_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile", 
                       mo_ohMPHP_i_cor_t3_ter == "2nd tertile" & mo_ohMPHP_i_cor_t3_ter_bis == "2nd tertile" ~ "2nd tertile", 
                       mo_ohMPHP_i_cor_t3_ter == "3rd tertile" & mo_ohMPHP_i_cor_t3_ter_bis == "3rd tertile" ~ "3rd tertile", 
                       
                       mo_ohMPHP_i_cor_t3_ter == "1st tertile" & mo_ohMPHP_i_cor_t3_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                       mo_ohMPHP_i_cor_t3_ter == "2nd tertile" & mo_ohMPHP_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                       
                       mo_ohMPHP_i_cor_t3_ter == "1st tertile" & mo_ohMPHP_i_cor_t3_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                       mo_ohMPHP_i_cor_t3_ter == "3rd tertile" & mo_ohMPHP_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                       
                       mo_ohMPHP_i_cor_t3_ter == "2nd tertile" & mo_ohMPHP_i_cor_t3_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                       mo_ohMPHP_i_cor_t3_ter == "3rd tertile" & mo_ohMPHP_i_cor_t3_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    
    DINCH = case_when(mo_DINCH_ms_i_cor_t3_ter == "1st tertile" & mo_DINCH_ms_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile", 
                      mo_DINCH_ms_i_cor_t3_ter == "2nd tertile" & mo_DINCH_ms_i_cor_t3_ter_bis == "2nd tertile" ~ "2nd tertile", 
                      mo_DINCH_ms_i_cor_t3_ter == "3rd tertile" & mo_DINCH_ms_i_cor_t3_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                      mo_DINCH_ms_i_cor_t3_ter == "1st tertile" & mo_DINCH_ms_i_cor_t3_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                      mo_DINCH_ms_i_cor_t3_ter == "2nd tertile" & mo_DINCH_ms_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                      mo_DINCH_ms_i_cor_t3_ter == "1st tertile" & mo_DINCH_ms_i_cor_t3_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                      mo_DINCH_ms_i_cor_t3_ter == "3rd tertile" & mo_DINCH_ms_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                      mo_DINCH_ms_i_cor_t3_ter == "2nd tertile" & mo_DINCH_ms_i_cor_t3_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                      mo_DINCH_ms_i_cor_t3_ter == "3rd tertile" & mo_DINCH_ms_i_cor_t3_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"))%>%
  select(ident, ident_bis, Bray_curtis_dissimilarity, 
         DEHP, MnBP, DiNP, MiBP, MBzP, MEP, ohMPHP, DINCH)



# on fait passer les données en long 
bdd_long_t3 <- bdd_long_t3 %>%
  pivot_longer(cols = -c("ident", "ident_bis", "Bray_curtis_dissimilarity"), 
               names_to = "Pollutant", 
               values_to = "Groups")


# duplication de 1st tertile, 2nd tertile, 3rd tertile
test <- 
  bdd_long_t3 %>% 
  filter(Groups %in% c("1st tertile", "2nd tertile", "3rd tertile")) %>%  
  filter(!Pollutant == "BPS") %>%
  mutate(Groups = fct_recode(Groups, 
                             " 1st tertile" = "1st tertile", 
                             " 2nd tertile" = "2nd tertile", 
                             " 3rd tertile" = "3rd tertile"))
bdd_long_t3 <- bind_rows(bdd_long_t3, test)
rm(test)

bdd_long_t3 <- bdd_long_t3 %>%
  mutate(
    Pollutant_rec = fct_recode(Pollutant, 
                               "ΣDEHP trim.3" = "DEHP",
                               "ΣDINCH trim.3" = "DINCH",
                               "ΣDiNP trim.3" = "DiNP",
                               "MBzP trim.3" = "MBzP",
                               "MEP trim.3" = "MEP",
                               "MiBP trim.3" = "MiBP",
                               "MnBP trim.3" = "MnBP",
                               "ohMPHP trim.3" = "ohMPHP"),
    Pollutant = fct_relevel(Pollutant,
                            "DEHP", "MnBP", "DiNP", "MiBP", "MBzP", "MEP", "ohMPHP", "DINCH"),
    Pollutant_rec = fct_relevel(Pollutant_rec,
                                "ΣDEHP trim.3", "MnBP trim.3", "ΣDiNP trim.3", "MiBP trim.3", 
                                "MBzP trim.3", "MEP trim.3", "ohMPHP trim.3", "ΣDINCH trim.3"),
    Groups = fct_relevel(Groups,
                         " 3rd tertile",                             # intra groupe forte exposition 
                         "2nd tertile - 3rd tertile",             # inter groupe moyenne et forte exposition
                         " 2nd tertile",                                  # intra groupe moyenne exposition 
                         
                         "3rd tertile",                                    # intra groupe forte exposition 
                         "1st tertile - 3rd tertile", # inter groupe faible et forte exposition 
                         " 1st tertile",
                         
                         "2nd tertile",                                  # intra groupe moyenne exposition 
                         "1st tertile - 2nd tertile",             # inter groupe faible et moyenne exposition 
                         "1st tertile"),                               # intra groupe faible exposition
    Groups_rec = fct_recode(Groups, 
                            "Medium vs High" = " 3rd tertile",
                            "Medium vs High" = "2nd tertile - 3rd tertile",
                            "Medium vs High" = " 2nd tertile",
                            "Low vs High" = "3rd tertile",
                            "Low vs High" = "1st tertile - 3rd tertile",
                            "Low vs High" = " 1st tertile",
                            "Low vs Medium" = "2nd tertile",
                            "Low vs Medium" = "1st tertile - 2nd tertile",
                            "Low vs Medium" = "1st tertile"))



#### boxplot version 1 ----
boxplot_t3_phthalates <- bdd_long_t3 %>%                            # intra groupe faible exposition
  ggplot() +
  aes(
    x = Bray_curtis_dissimilarity,
    y = Groups,
    fill = Groups
  ) +
  geom_boxplot() +  
  scale_fill_manual(
    values = c(`1st tertile` = "#FFE0DE",   # les catégories faibles expositions (intragroupes)
               ` 1st tertile` = "#FFE0DE", 
               
               `1st tertile - 2nd tertile` = "gray70",  # les variances intergroupes
               `1st tertile - 3rd tertile` = "gray70",
               `2nd tertile - 3rd tertile` = "gray70",
               
               `2nd tertile` = "#FF8F87",   # les catégorires moyennes expositions (intragroupes)
               ` 2nd tertile` = "#FF8F87",   
               
               `3rd tertile` = "#FF4034",   
               ` 3rd tertile` = "#FF4034")
  ) +
  labs(x = "Bray Curtis dissimilarity", 
       y = "") +
  theme_lucid() +
  theme(legend.position = "none", 
        axis.title = element_text(face = "bold")) +
  facet_wrap(~Pollutant_rec, scales = "free", ncol = 2)




### 1 year ----
# on réduit la matrice de dissimilarité pour ne pas avoir les données en double
metric_bray_curtis_red <- metric_bray_curtis %>% column_to_rownames(var = "ident")
metric_bray_curtis_red[lower.tri(metric_bray_curtis_red)]<- NA
metric_bray_curtis_red <- metric_bray_curtis_red %>% rownames_to_column(var = "ident")

# on merge les données d'expo t3 aux données de dissimilarité
bdd_long_Y1 <- bdd_final_Y1 %>% select(ident, all_of(explanatory_vars_Y1))
bdd_long_Y1 <- full_join(bdd_long_Y1, metric_bray_curtis_red, by = "ident")

# on fait passer la base de donnnées en long
bdd_long_Y1[,explanatory_vars_Y1] <- lapply(bdd_long_Y1[,explanatory_vars_Y1], as.character)
bdd_long_Y1 <- bdd_long_Y1 %>%
  pivot_longer(cols = c(-"ident", -all_of(explanatory_vars_Y1)), 
               names_to = "ident_bis", 
               values_to = "Bray_curtis_dissimilarity") %>%
  select(ident, ident_bis, Bray_curtis_dissimilarity, everything()) %>%
  filter(ident != ident_bis) %>%
  filter(!is.na(Bray_curtis_dissimilarity)) 

# on créé une base de données bis pour savoir à quelle groupes d'expos sont comparés les paires 
bdd_long_Y1_bis <- bdd_final_Y1 %>%
  select(ident, 
         all_of(explanatory_vars_Y1)) %>%
  rename(ident_bis = ident, 
         ch_DEHP_ms_i_cor_Y1_ter_bis = ch_DEHP_ms_i_cor_Y1_ter, 
         ch_MnBP_i_cor_Y1_ter_bis = ch_MnBP_i_cor_Y1_ter, 
         ch_DiNP_ms_i_cor_Y1_ter_bis = ch_DiNP_ms_i_cor_Y1_ter, 
         ch_MiBP_i_cor_Y1_ter_bis = ch_MiBP_i_cor_Y1_ter, 
         ch_MBzP_i_cor_Y1_ter_bis = ch_MBzP_i_cor_Y1_ter, 
         ch_MEP_i_cor_Y1_ter_bis = ch_MEP_i_cor_Y1_ter, 
         ch_ohMPHP_i_cor_Y1_ter_bis = ch_ohMPHP_i_cor_Y1_ter, 
         ch_DINCH_ms_i_cor_Y1_ter_bis = ch_DINCH_ms_i_cor_Y1_ter)

# on merge les données puis on supprime la base de données bis 
bdd_long_Y1 <- full_join(bdd_long_Y1, bdd_long_Y1_bis, by = "ident_bis") 
bdd_long_Y1 <- bdd_long_Y1 %>%
  filter(!ident %in% c(23994, 25166, 26766, 14668, 26923)) %>%
  filter(!ident_bis %in% c(23994, 25166, 26766, 14668, 26923)) %>%
  filter(!is.na(ident))
rm(bdd_long_Y1_bis)

# on créé une variable qui indique de quel intergroupe il s'agit
bdd_long_Y1 <- bdd_long_Y1 %>%
  mutate(
    DEHP = case_when(ch_DEHP_ms_i_cor_Y1_ter == "1st tertile" & ch_DEHP_ms_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile", 
                     ch_DEHP_ms_i_cor_Y1_ter == "2nd tertile" & ch_DEHP_ms_i_cor_Y1_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     ch_DEHP_ms_i_cor_Y1_ter == "3rd tertile" & ch_DEHP_ms_i_cor_Y1_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     ch_DEHP_ms_i_cor_Y1_ter == "1st tertile" & ch_DEHP_ms_i_cor_Y1_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     ch_DEHP_ms_i_cor_Y1_ter == "2nd tertile" & ch_DEHP_ms_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     ch_DEHP_ms_i_cor_Y1_ter == "1st tertile" & ch_DEHP_ms_i_cor_Y1_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     ch_DEHP_ms_i_cor_Y1_ter == "3rd tertile" & ch_DEHP_ms_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     ch_DEHP_ms_i_cor_Y1_ter == "2nd tertile" & ch_DEHP_ms_i_cor_Y1_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     ch_DEHP_ms_i_cor_Y1_ter == "3rd tertile" & ch_DEHP_ms_i_cor_Y1_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    MnBP = case_when(ch_MnBP_i_cor_Y1_ter == "1st tertile" & ch_MnBP_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile", 
                     ch_MnBP_i_cor_Y1_ter == "2nd tertile" & ch_MnBP_i_cor_Y1_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     ch_MnBP_i_cor_Y1_ter == "3rd tertile" & ch_MnBP_i_cor_Y1_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     ch_MnBP_i_cor_Y1_ter == "1st tertile" & ch_MnBP_i_cor_Y1_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     ch_MnBP_i_cor_Y1_ter == "2nd tertile" & ch_MnBP_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     ch_MnBP_i_cor_Y1_ter == "1st tertile" & ch_MnBP_i_cor_Y1_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     ch_MnBP_i_cor_Y1_ter == "3rd tertile" & ch_MnBP_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     ch_MnBP_i_cor_Y1_ter == "2nd tertile" & ch_MnBP_i_cor_Y1_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     ch_MnBP_i_cor_Y1_ter == "3rd tertile" & ch_MnBP_i_cor_Y1_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    DiNP = case_when(ch_DiNP_ms_i_cor_Y1_ter == "1st tertile" & ch_DiNP_ms_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile", 
                     ch_DiNP_ms_i_cor_Y1_ter == "2nd tertile" & ch_DiNP_ms_i_cor_Y1_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     ch_DiNP_ms_i_cor_Y1_ter == "3rd tertile" & ch_DiNP_ms_i_cor_Y1_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     ch_DiNP_ms_i_cor_Y1_ter == "1st tertile" & ch_DiNP_ms_i_cor_Y1_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     ch_DiNP_ms_i_cor_Y1_ter == "2nd tertile" & ch_DiNP_ms_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     ch_DiNP_ms_i_cor_Y1_ter == "1st tertile" & ch_DiNP_ms_i_cor_Y1_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     ch_DiNP_ms_i_cor_Y1_ter == "3rd tertile" & ch_DiNP_ms_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     ch_DiNP_ms_i_cor_Y1_ter == "2nd tertile" & ch_DiNP_ms_i_cor_Y1_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     ch_DiNP_ms_i_cor_Y1_ter == "3rd tertile" & ch_DiNP_ms_i_cor_Y1_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    MiBP = case_when(ch_MiBP_i_cor_Y1_ter == "1st tertile" & ch_MiBP_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile", 
                     ch_MiBP_i_cor_Y1_ter == "2nd tertile" & ch_MiBP_i_cor_Y1_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     ch_MiBP_i_cor_Y1_ter == "3rd tertile" & ch_MiBP_i_cor_Y1_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     ch_MiBP_i_cor_Y1_ter == "1st tertile" & ch_MiBP_i_cor_Y1_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     ch_MiBP_i_cor_Y1_ter == "2nd tertile" & ch_MiBP_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     ch_MiBP_i_cor_Y1_ter == "1st tertile" & ch_MiBP_i_cor_Y1_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     ch_MiBP_i_cor_Y1_ter == "3rd tertile" & ch_MiBP_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     ch_MiBP_i_cor_Y1_ter == "2nd tertile" & ch_MiBP_i_cor_Y1_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     ch_MiBP_i_cor_Y1_ter == "3rd tertile" & ch_MiBP_i_cor_Y1_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    MBzP = case_when(ch_MBzP_i_cor_Y1_ter == "1st tertile" & ch_MBzP_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile", 
                     ch_MBzP_i_cor_Y1_ter == "2nd tertile" & ch_MBzP_i_cor_Y1_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     ch_MBzP_i_cor_Y1_ter == "3rd tertile" & ch_MBzP_i_cor_Y1_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     ch_MBzP_i_cor_Y1_ter == "1st tertile" & ch_MBzP_i_cor_Y1_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     ch_MBzP_i_cor_Y1_ter == "2nd tertile" & ch_MBzP_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     ch_MBzP_i_cor_Y1_ter == "1st tertile" & ch_MBzP_i_cor_Y1_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     ch_MBzP_i_cor_Y1_ter == "3rd tertile" & ch_MBzP_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     ch_MBzP_i_cor_Y1_ter == "2nd tertile" & ch_MBzP_i_cor_Y1_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     ch_MBzP_i_cor_Y1_ter == "3rd tertile" & ch_MBzP_i_cor_Y1_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    MEP = case_when(ch_MEP_i_cor_Y1_ter == "1st tertile" & ch_MEP_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile", 
                    ch_MEP_i_cor_Y1_ter == "2nd tertile" & ch_MEP_i_cor_Y1_ter_bis == "2nd tertile" ~ "2nd tertile", 
                    ch_MEP_i_cor_Y1_ter == "3rd tertile" & ch_MEP_i_cor_Y1_ter_bis == "3rd tertile" ~ "3rd tertile", 
                    
                    ch_MEP_i_cor_Y1_ter == "1st tertile" & ch_MEP_i_cor_Y1_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                    ch_MEP_i_cor_Y1_ter == "2nd tertile" & ch_MEP_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                    
                    ch_MEP_i_cor_Y1_ter == "1st tertile" & ch_MEP_i_cor_Y1_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                    ch_MEP_i_cor_Y1_ter == "3rd tertile" & ch_MEP_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                    
                    ch_MEP_i_cor_Y1_ter == "2nd tertile" & ch_MEP_i_cor_Y1_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                    ch_MEP_i_cor_Y1_ter == "3rd tertile" & ch_MEP_i_cor_Y1_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    ohMPHP = case_when(ch_ohMPHP_i_cor_Y1_ter == "1st tertile" & ch_ohMPHP_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile", 
                       ch_ohMPHP_i_cor_Y1_ter == "2nd tertile" & ch_ohMPHP_i_cor_Y1_ter_bis == "2nd tertile" ~ "2nd tertile", 
                       ch_ohMPHP_i_cor_Y1_ter == "3rd tertile" & ch_ohMPHP_i_cor_Y1_ter_bis == "3rd tertile" ~ "3rd tertile", 
                       
                       ch_ohMPHP_i_cor_Y1_ter == "1st tertile" & ch_ohMPHP_i_cor_Y1_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                       ch_ohMPHP_i_cor_Y1_ter == "2nd tertile" & ch_ohMPHP_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                       
                       ch_ohMPHP_i_cor_Y1_ter == "1st tertile" & ch_ohMPHP_i_cor_Y1_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                       ch_ohMPHP_i_cor_Y1_ter == "3rd tertile" & ch_ohMPHP_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                       
                       ch_ohMPHP_i_cor_Y1_ter == "2nd tertile" & ch_ohMPHP_i_cor_Y1_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                       ch_ohMPHP_i_cor_Y1_ter == "3rd tertile" & ch_ohMPHP_i_cor_Y1_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    DINCH = case_when(ch_DINCH_ms_i_cor_Y1_ter == "1st tertile" & ch_DINCH_ms_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile", 
                     ch_DINCH_ms_i_cor_Y1_ter == "2nd tertile" & ch_DINCH_ms_i_cor_Y1_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     ch_DINCH_ms_i_cor_Y1_ter == "3rd tertile" & ch_DINCH_ms_i_cor_Y1_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     ch_DINCH_ms_i_cor_Y1_ter == "1st tertile" & ch_DINCH_ms_i_cor_Y1_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     ch_DINCH_ms_i_cor_Y1_ter == "2nd tertile" & ch_DINCH_ms_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     ch_DINCH_ms_i_cor_Y1_ter == "1st tertile" & ch_DINCH_ms_i_cor_Y1_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     ch_DINCH_ms_i_cor_Y1_ter == "3rd tertile" & ch_DINCH_ms_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     ch_DINCH_ms_i_cor_Y1_ter == "2nd tertile" & ch_DINCH_ms_i_cor_Y1_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     ch_DINCH_ms_i_cor_Y1_ter == "3rd tertile" & ch_DINCH_ms_i_cor_Y1_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"))%>%
  select(ident, ident_bis, Bray_curtis_dissimilarity, 
         DEHP, MnBP, DiNP, MiBP, MBzP, MEP, ohMPHP, DINCH)



# on fait passer les données en long 
bdd_long_Y1 <- bdd_long_Y1 %>%
  pivot_longer(cols = -c("ident", "ident_bis", "Bray_curtis_dissimilarity"), 
               names_to = "Pollutant", 
               values_to = "Groups")


# duplication de 1st tertile, 2nd tertile, 3rd tertile
test <- 
  bdd_long_Y1 %>% 
  filter(Groups %in% c("1st tertile", "2nd tertile", "3rd tertile")) %>%  
  mutate(Groups = fct_recode(Groups, 
                             " 1st tertile" = "1st tertile", 
                             " 2nd tertile" = "2nd tertile", 
                             " 3rd tertile" = "3rd tertile"))
bdd_long_Y1 <- bind_rows(bdd_long_Y1, test)
rm(test)

bdd_long_Y1 <- bdd_long_Y1 %>%
    mutate(
      Pollutant_rec = fct_recode(Pollutant, 
                                 "ΣDEHP 12 months" = "DEHP",
                                 "ΣDINCH 12 months" = "DINCH",
                                 "ΣDiNP 12 months" = "DiNP",
                                 "MBzP 12 months" = "MBzP",
                                 "MEP 12 months" = "MEP",
                                 "MiBP 12 months" = "MiBP",
                                 "MnBP 12 months" = "MnBP",
                                 "ohMPHP 12 months" = "ohMPHP"),
      Pollutant = fct_relevel(Pollutant,
                              "DEHP", "MnBP", "DiNP", "MiBP", "MBzP", "MEP", "ohMPHP", "DINCH"),
      Pollutant_rec = fct_relevel(Pollutant_rec,
                                  "ΣDEHP 12 months", "MnBP 12 months", "ΣDiNP 12 months", "MiBP 12 months", 
                                  "MBzP 12 months", "MEP 12 months", "ohMPHP 12 months", "ΣDINCH 12 months"),
    Groups = fct_relevel(Groups,
                         " 3rd tertile",                           # intra groupe forte exposition 
                         "2nd tertile - 3rd tertile",             # inter groupe moyenne et forte exposition
                         " 2nd tertile",                                   # intra groupe moyenne exposition 
                         
                         "3rd tertile",                                    # intra groupe forte exposition 
                         "1st tertile - 3rd tertile", # inter groupe faible et forte exposition 
                         " 1st tertile",
                         
                         "2nd tertile",                                  # intra groupe moyenne exposition 
                         "1st tertile - 2nd tertile",             # inter groupe faible et moyenne exposition 
                         "1st tertile"))                               # intra groupe faible exposition



#### boxplot version 1 ----
boxplot_Y1_phthalates <- bdd_long_Y1 %>%                            # intra groupe faible exposition
  ggplot() +
  aes(
    x = Bray_curtis_dissimilarity,
    y = Groups,
    fill = Groups
  ) +
  geom_boxplot() +  
  scale_fill_manual(
    values = c(`1st tertile` = "#FFE0DE",   # les catégories faibles expositions (intragroupes)
               ` 1st tertile` = "#FFE0DE", 
               
               `1st tertile - 2nd tertile` = "gray70",  # les variances intergroupes
               `1st tertile - 3rd tertile` = "gray70",
               `2nd tertile - 3rd tertile` = "gray70",
               
               
               `2nd tertile` = "#FF8F87",   # les catégorires moyennes expositions (intragroupes)
               ` 2nd tertile` = "#FF8F87",   
               
               `3rd tertile` = "#FF4034",  
               ` 3rd tertile` = "#FF4034")
  ) +
  labs(x = "Bray Curtis dissimilarity", 
       y = "") +
  theme_lucid() +
  theme(legend.position = "none", 
        axis.title = element_text(face = "bold")) +
  facet_wrap(~Pollutant_rec, scales = "free", ncol = 2)

ggsave(plot = boxplot_t2_phthalates, 
       filename = "4_output/betadiv/boxplot_t2_phthalates.tiff", 
       device = "tiff", 
       units = "cm", 
       width = 25, 
       height = 25)

ggsave(plot = boxplot_t3_phthalates, 
       filename = "4_output/betadiv/boxplot_t3_phthalates.tiff", 
       device = "tiff", 
       units = "cm", 
       width = 25, 
       height = 25)

ggsave(plot = boxplot_Y1_phthalates, 
       filename = "4_output/betadiv/boxplot_Y1_phthalates.tiff", 
       device = "tiff", 
       units = "cm", 
       width = 25, 
       height = 25)


# Test 1 seule fois vérif ----
set.seed(1996)
results <- 
  adonis2(
    all_dist_t2 ~ mo_DEHP_ms_i_cor_t2_ter,
    data = bdd_final_t2,
    permutations = 999)

results <- 
  adonis2(
    all_dist_t2 ~ 
      mo_DEHP_ms_i_cor_t2_ter +
      ch_feces_RUN_Y1 +
      ch_feces_age_w_Y1_i +
      po_delmod +
      ch_food_intro_Y1_3cat_i +
      ch_antibio_Y1_2cat_i +
      mo_par_2cat +
      mo_pets_i +
      ch_sex +
      mo_tob_gr_anyt_yn_n2_i +
      Mo_ETS_anyT_yn1_opt_i + 
      ch_ETS_12m_opt36m +
      # mo_interpreg_3cat +
      mo_dipl_3cat_i +
      po_w_kg_3cat +
      po_he_3cat_i +
      ch_w_Y1_3cat_i +
      ch_he_Y1_3cat_i +
      po_gd +
      mo_age +
      mo_bmi_bepr_3cat_i +
      bf_duration_till48w_4cat_i,
    data = bdd_final_t2,
    permutations = 999
  )



