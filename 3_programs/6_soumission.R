## Aline Davias
## Analyses microbiot et neurodeveloppement 
## 20/02/2024


# Packages loading ----
library(tidyverse)
library(phyloseq)
library(expss)
library(gtsummary)
library(labelled)
library(questionr)
library(sjlabelled)
library(openxlsx)
library(corrplot)
library(see)
library(psych)
library(compositions)
library(writexl)
library(performance)
library(see)
library(broom)
library(broom.helpers)
library(broom.mixed)
library(mice)
library(gt)
library(GGally)
theme_gtsummary_compact()

# Data and functions loading ----
load("1_intermediate_data/3_data_imputation_AD_gumme.RData")
source("3_programs/4_functions_AD_gumme.R")
rm(table_cor, table_cor_sg, model_Y3_1, model_Y3_2, heatmap_cor, heatmap_cor_pairwise)
corres <- 
  taxa_table_Y1 %>% 
  select(Phyla_corres = ch_feces_phylum_ASVbased_Y1, 
         Class_corres = ch_feces_class_ASVbased_Y1,
         Order_corres = ch_feces_order_ASVbased_Y1, 
         Family_corres = ch_feces_family_ASVbased_Y1,
         Exposure = ch_feces_genus_ASVbased_Y1) %>%
  filter(Exposure %in% genera_linear) %>%
  distinct(Exposure, .keep_all = TRUE)

bdd_final_imp_1 <- bdd_final_imp_1 %>%
  mutate(ch_feces_rel_p1_Y1_10 = ch_feces_rel_p1_Y1/10, 
         ch_feces_rel_p2_Y1_10 = ch_feces_rel_p2_Y1/10, 
         ch_feces_rel_p3_Y1_10 = ch_feces_rel_p3_Y1/10, 
         ch_feces_rel_p4_Y1_10 = ch_feces_rel_p4_Y1/10)

# Vectors ----
covariates <- covar_vec_model_final
explanatory <- bdd_final_imp_1 %>% 
  select(all_of(microbiote_vec), 
         ch_feces_rel_p1_Y1_10, 
         ch_feces_rel_p2_Y1_10, 
         ch_feces_rel_p3_Y1_10, 
         ch_feces_rel_p4_Y1_10) %>% 
  select(-ch_feces_SpecRich_cmin_ASV_Y1, 
         -ch_feces_Shannon_cmin_ASV_Y1,
         -ch_feces_SpecRich_10000_ASV_Y1, 
         -ch_feces_Shannon_10000_ASV_Y1, 
         -ch_feces_rel_p1_Y1,
         -ch_feces_rel_p2_Y1,
         -ch_feces_rel_p3_Y1,
         -ch_feces_rel_p4_Y1) %>%
  colnames()
outcomes <- bdd %>%
  select(all_of(neuro_vec)) %>%
  select(-ch_socawar_y3,
         -ch_soccog_y3,
         -ch_soccom_y3, 
         -ch_socmot_y3, 
         -ch_RRB_y3) %>%
  colnames()
rm(neuro_vec, microbiote_vec, covar_a_imputer, covar_a_tester, covar_vec_model_final)

genera_linear_complet <- bdd %>%
  select(contains("ch_feces_rel_g")) %>%
  filter(!is.na(ch_feces_rel_g1_Y1)) %>%
  select_if(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.3) %>%
  colnames()

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

# Data description ----
## Covariates ----
descrip_covar <- 
  bdd_final_imp_1 %>% 
  select(all_of(covariates)) %>% 
  tbl_summary(type = list(all_categorical() ~ "categorical"), 
              label = list(ch_care_main_12m_opt2_2c ~ "Main mode of child care at 12 months")) %>%
  bold_labels()

descrip_covar

## Included genera (n=46) ----
### Tables ----
bdd %>% 
  select(all_of(genera_linear_complet)) %>%
  na.omit()%>%
  View()

descrip_genera_linear <- tbl_merge(
  tbls = 
    list(
      tbl_1 = 
        bdd %>% 
        select(all_of(genera_linear_complet)) %>%
        na.omit()%>%
        set_label(genera_linear) %>%
        tbl_summary(
          type = list(everything() ~ "continuous"), 
          statistic = list(everything() ~ "{median} ({p25}, {p75})"), 
          digits = list(all_continuous() ~ c(2, 1, 1))), 
      tbl_2 = 
        bdd_final_imp_1 %>% 
        select(all_of(genera_linear)) %>%
        set_label(genera_linear) %>%
        tbl_summary(
          type = list(everything() ~ "continuous"), 
          statistic = list(everything() ~ "{median} ({p25}, {p75})"), 
          digits = list(all_continuous() ~ c(2, 1, 1))), 
      tbl_3 = 
        bdd %>% 
        select(all_of(genera_linear_complet)) %>%
        na.omit()%>%
        mutate_all(~ ifelse(.>0, "Yes", "No")) %>%
        set_label(genera_linear) %>%
        tbl_summary(
          type = list(everything() ~ "categorical"))), 
  tab_spanner = c("**Continuous**", "**Log transformed**", "**Categorical (Y/N)**"))
descrip_genera_linear

### Density plots ----
#### Distribution avant imputation et transformation 
densityplot(data = bdd, vars = genera_linear_complet[1:23], ncol = 5)    
densityplot(data = bdd, vars = genera_linear_complet[24:46], ncol = 5)

#### Distribution après imputation et transformation 
densityplot(data = bdd_final_imp_1, vars = genera_linear[1:23], ncol = 5)    
densityplot(data = bdd_final_imp_1, vars = genera_linear[24:46], ncol = 5)


# Calcul du nombre de tests effectifs pour la correction pour comp multiple (FWER) ----
test_expo_taxa <- bdd_final_imp_1 %>%
  filter(!is.na(ch_feces_rel_p1_Y1)) %>%
  select(all_of(explanatory)) %>%
  select(-"ch_feces_SpecRich_5000_ASV_Y1", -"ch_feces_Shannon_5000_ASV_Y1") 
test_expo_taxa <- lapply(test_expo_taxa, function(x) { var_lab(x) <- NULL; return(x) })
test_expo_taxa <- as.data.frame(test_expo_taxa)
cor_mixed_table_expo_taxa <- cor_mixed(data = test_expo_taxa)
results_alpha_corrected_expo_taxa <- alpha_corrected(data = test_expo_taxa, alpha = 0.05)
results_M0_corrected_expo_taxa <- M0_corrected(data = test_expo_taxa, alpha = 0.05)
results_M0_corrected_expo_taxa

test_expo_outcomes <- bdd_final_imp_1 %>%
  select(all_of(outcomes))
test_expo_outcomes <- lapply(test_expo_outcomes, function(x) { var_lab(x) <- NULL; return(x) })
test_expo_outcomes <- as.data.frame(test_expo_outcomes)
cor_mixed_table_expo_outcomes <- cor_mixed(data = test_expo_outcomes)
results_alpha_corrected_expo_outcomes <- alpha_corrected(data = test_expo_outcomes, alpha = 0.05)
results_M0_corrected_expo_outcomes <- M0_corrected(data = test_expo_outcomes, alpha = 0.05)
results_M0_corrected_expo_outcomes

## Pour les analyses de diversité, le seuil corrigé est : 0.05/(3*7)
## Pour les analyses de taxonomie, le seuil corrigé est : 0.05/(33*7)
rm(test_expo_taxa, test_expo_outcomes, 
   cor_mixed_table_expo_taxa, cor_mixed_table_expo_outcomes, 
   results_alpha_corrected_expo_taxa, results_alpha_corrected_expo_outcomes,
   results_M0_corrected_expo_taxa, results_M0_corrected_expo_outcomes, 
   cor_mixed, alpha_corrected, M0_corrected)

# Running linear regressions ----
## adapation des covariables selon l'outcome
covariates_CBCL <- c("ch_age_CBCL_Y2", covariates)                    
covariates_IQ <- c("ch_age_IQ_Y3", covariates)
covariates_SRS_BRIEF <- c("ch_age_SRS_BRIEFP_Y3", covariates)

covariates_map <- list(
  CBCL = covariates_CBCL,
  SRS_BRIEF = covariates_SRS_BRIEF,
  IQ = covariates_IQ
)

## Création d'une liste pour stocker les tableaux par outcome
tbls_by_outcome_multi <- vector("list", length(outcomes))
names(tbls_by_outcome_multi) <- outcomes

# Boucle principale
for (outcome in outcomes) {
  # Sélection des covariables appropriées
  if (outcome %in% c("ch_cbclintscore_y2", "ch_cbclextscore_y2")) {
    covariates <- covariates_map$CBCL
  } else if (outcome %in% c("ch_SRStotal_y3", "ch_briefpinhibit_y3", "ch_briefpshift_y3", "ch_briefpemocontrol_y3", "ch_briefpworkmemo_y3", "ch_briefpplan_y3")) {
    covariates <- covariates_map$SRS_BRIEF
  } else if (outcome %in% c("ch_verbal_comprehension_IQ_Y3", "ch_visuospatiale_IQ_Y3", "ch_work_memory_IQ_Y3", "ch_total_IQ_Y3")) {
    covariates <- covariates_map$IQ
  }
  
  tbls_for_outcome_multi <- vector("list", length(explanatory))
  names(tbls_for_outcome_multi) <- explanatory
  
  for (exposure in explanatory) {                                               # running linear regression
    terms <- c(exposure, covariates)
    formula <- reformulate(terms, response = outcome)
    model <- lm(formula, data = bdd_final_imp_1)
    # model <- with(data = bdd_final_imp, 
    #               exp = lm(formula))
    
    tbl <-                                                                      
      tbl_regression(
        model, 
        include = exposure,
        estimate_fun = scales::label_number(accuracy = .01, decimal.mark = "."),
        pvalue_fun = custom_pvalue_fun,
        exponentiate = FALSE) %>%
      bold_p() %>%
      bold_labels() %>%
      add_global_p(include = exposure, singular.ok = TRUE, keep = TRUE)
    
    tbls_for_outcome_multi[[exposure]] <- tbl
  }
  tbls_by_outcome_multi[[outcome]] <- tbls_for_outcome_multi
}

table_multi <- tibble()                                                         # Initialisation d'un tibble vide pour stocker les résultats finaux
for (i in seq_along(tbls_by_outcome_multi)) {                                   # Nom de l'outcome pour cette itération
  outcome_name <- names(tbls_by_outcome_multi)[i]
  
  for (j in seq_along(tbls_by_outcome_multi[[i]])) {                            # Itération sur chaque tbl_regression dans la liste courante
    exposure_name <- names(tbls_by_outcome_multi[[i]])[j]                       # Nom de la variable d'exposition pour cette itération
    tbl_data <- tbls_by_outcome_multi[[i]][[j]] %>%                             # Extraction des données du tableau tbl_regression
      as_tibble() %>%
      mutate(Outcome = outcome_name, Exposure = exposure_name)
    
    table_multi <- bind_rows(table_multi, tbl_data)                             # Ajout des données extraites au tibble final
  }
}
rm(terms, formula, model,
   tbl, tbl_data,
   tbls_for_outcome_multi,
   exposure, exposure_name, i, j, outcome, outcome_name)
covariates <- c("po_w_kg_3cat",                                                 # redéfinition de covariates sinon il garde une variable age 
                "po_he_3cat",
                "mo_dipl_2cat",
                "mo_age",
                "mo_bmi_bepr_3cat",
                "ch_sex",
                "mo_par_2cat",
                "ch_bf_duration_till48w_4cat",
                "po_gd",
                "po_delmod",
                "ch_food_intro_Y1_3cat",
                "mo_pets",
                "ch_antibio_Y1_2cat",
                "home_total_y3",              
                "mo_hadtotscore_grt3_imp",
                "mo_tob_gr_anyt_yn_n2",
                "ch_tabacco_passive_up_to_Y1",
                "ch_care_main_12m_opt2_2c")

table_multi <-                                                                         # Ajout des correspondances taxonomiques 
  left_join(table_multi, corres, by = "Exposure") %>%
  select(Phyla_corres, Class_corres, Order_corres, Family_corres, everything())
rm(corres)

table_multi <- table_multi %>%
  select(Outcome, 
         Phyla_corres, Class_corres, Order_corres, Family_corres,
         `Gut microbiota parameters` = Exposure, 
         Beta = "**Beta**", 
         "95% CI" = "**95% CI**",
         "p-value" = "**p-value**", 
         "Characteristic" = "**Characteristic**") %>%
  mutate(
    `Gut microbiota parameters` = 
      fct_recode(`Gut microbiota parameters`, 
                 "Firmicutes" = "ch_feces_rel_p1_Y1_10",
                 "Actinobacteria" = "ch_feces_rel_p2_Y1_10",
                 "Bacteroidetes" = "ch_feces_rel_p3_Y1_10",
                 "Proteobacteria" = "ch_feces_rel_p4_Y1_10",
                 "Shannon diversity" = "ch_feces_Shannon_5000_ASV_Y1",
                 "Specific richness" = "ch_feces_SpecRich_5000_ASV_Y1",
                 "Clostridium IV" = "Clostridium_IV",
                  "Clostridium sensu stricto" = "Clostridium_sensu_stricto",
                  "Clostridium XlVa" = "Clostridium_XlVa",
                  "Clostridium XVIII" = "Clostridium_XVIII",
                  "Erysipelotrichaceae incertae sedis" = "Erysipelotrichaceae_incertae_sedis",
                  "Escherichia and Shigella" = "Escherichia_Shigella",
                  "Lachnospiracea incertae sedis" = "Lachnospiracea_incertae_sedis",
                  "Ruminococcus 2" = "Ruminococcus2",
                  "Saccharibacteria genera incertae sedis" = "Saccharibacteria_genera_incertae_sedis"),
    `Gut microbiota parameters` = 
      fct_relevel(`Gut microbiota parameters`, 
                  "Saccharibacteria genera incertae sedis",
                  "Peptoniphilus", "Granulicatella", "Anaerotruncus", "Lactococcus", 
                  "Terrisporobacter", "Oscillibacter", "Haemophilus",
                  "Coprococcus", "Erysipelotrichaceae incertae sedis", "Butyricicoccus", 
                  "Dialister",  "Subdoligranulum", "Intestinibacter", "Klebsiella",
                  "Eisenbergiella", "Hungatella","Dorea", "Eggerthella",
                  "Romboutsia","Clostridium IV","Ruminococcus 2","Flavonifractor", "Alistipes",
                  "Collinsella","Parabacteroides","Fusicatenibacter", "Veillonella",
                  "Roseburia", "Enterobacter","Cellulosibacter",  "Enterococcus",
                  "Clostridium sensu stricto", "Clostridium XVIII", "Ruminococcus", 
                   "Gemmiger", "Anaerostipes", "Lachnospiracea incertae sedis", 
                  "Clostridium XlVa",  "Streptococcus", "Faecalibacterium", "Akkermansia",
                  "Escherichia and Shigella", "Blautia", "Bacteroides","Bifidobacterium", 
                  "Proteobacteria", "Bacteroidetes","Actinobacteria", "Firmicutes",
                  "Shannon diversity", "Specific richness"),
    Outcome = 
      fct_recode(Outcome, 
                 "Emotional control BRIEF-P score at 3 years" = "ch_briefpemocontrol_y3",
                  "Inhibition BRIEF-P score at 3 years" = "ch_briefpinhibit_y3",
                  "Plan and organization BRIEF-P score at 3 years" = "ch_briefpplan_y3",
                  "Shift BRIEF-P score at 3 years" = "ch_briefpshift_y3",
                  "Working memory BRIEF-P score at 3 years" = "ch_briefpworkmemo_y3",
                  "Externalizing CBCL score at 2 years" = "ch_cbclextscore_y2",
                  "Internalizing CBCL score at 2 years" = "ch_cbclintscore_y2",
                  "Total SRS score at 3 years" = "ch_SRStotal_y3",
                  "Total IQ score at 3 years" = "ch_total_IQ_Y3",
                  "Verbal comprehension IQ score at 3 years" = "ch_verbal_comprehension_IQ_Y3",
                  "Visuospatiale IQ score at 3 years" = "ch_visuospatiale_IQ_Y3",
                  "Work memory IQ score at 3 years" = "ch_work_memory_IQ_Y3"),
    Outcome = 
      fct_relevel(Outcome,
                  "Internalizing CBCL score at 2 years", "Externalizing CBCL score at 2 years",
                  "Emotional control BRIEF-P score at 3 years", "Inhibition BRIEF-P score at 3 years",
                  "Plan and organization BRIEF-P score at 3 years", "Shift BRIEF-P score at 3 years",
                  "Working memory BRIEF-P score at 3 years", "Total SRS score at 3 years",
                  "Verbal comprehension IQ score at 3 years", "Visuospatiale IQ score at 3 years",
                  "Work memory IQ score at 3 years", "Total IQ score at 3 years"),
    `p-value` = gsub("__", "", `p-value`),
    `p-value` = as.numeric(`p-value`),
    `q-value` = case_when(
      `Gut microbiota parameters` == "Specific richness" ~ `p-value`/(3*7),     # correction spéciale analyses de diversité 
      `Gut microbiota parameters` == "Shannon diversity" ~ `p-value`/(3*7),     # correction spéciale analyses de diversité 
      .default = `p-value`/(33*7)),                                             # correction spéciale analyses de taxonomie 
    p_value_shape = ifelse(`p-value`<0.05, "p-value<0.05", "p-value≥0.05"),
    q_value_shape = ifelse(`q-value`<0.05, "q-value<0.05", "q-value≥0.05"), 
    sens_beta = ifelse(Beta < 0, "Beta<0", "Beta≥0"), 
    sens_beta = fct_relevel(sens_beta, "Beta≥0", "Beta<0"))  %>% 
  separate(col = "95% CI", into = c("lower_CI", "upper_CI"), sep = ",", remove = FALSE) %>%
  mutate(
    lower_CI = as.numeric(lower_CI),
    upper_CI = as.numeric(upper_CI)
  ) %>%
  select(
    Phyla_corres, Class_corres, Order_corres, Family_corres,
    Outcome, 
    `Gut microbiota parameters`,
    Beta, sens_beta, 
    "95% CI", lower_CI, upper_CI, 
    "p-value", p_value_shape, 
    "q-value", q_value_shape)


# Tables ----
## Table 1: Covariables description ----
table_1 <- bdd %>% 
  select(all_of(covariates), 
         all_of(covar_age), 
         ch_feces_rel_p1_Y1) %>%
  mutate(
    statu = ifelse(!is.na(ch_feces_rel_p1_Y1), "included", "excluded"), 
    statu = fct_relevel(statu,  "included", "excluded")) %>%
  select(-ch_feces_rel_p1_Y1) %>%
  tbl_summary(
    type = list(all_categorical() ~ "categorical"), 
    label = list(ch_care_main_12m_opt2_2c ~ "Main mode of child care at 12 months"), 
    by="statu", 
    missing = "no", 
    digits = list(mo_age ~ 0, 
                  po_gd ~ 0, 
                  home_total_y3 ~ 0, 
                  mo_hadtotscore_grt3_imp ~ 0, 
                  ch_age_CBCL_Y2 ~ 0, 
                  ch_age_IQ_Y3 ~ 0, 
                  ch_age_SRS_BRIEFP_Y3 ~ 0)) %>%
  add_p(pvalue_fun = function(x) formatC(x, format = "f", digits = 3)) %>%
  bold_labels() %>%
  add_n()

## Table 2: Associations alphadiv and neurodev ----
table_2_long <-
  tbl_stack(tbls =
              lapply(tbls_by_outcome_multi, function(tbl_list) {
                tbl_merge(
                  list(
                    tbl_list[[1]] %>% add_n(),   # correspond à ch_feces_SpecRich_5000_ASV_Y1
                    tbl_list[[2]] %>% add_n()),  # correspond àch_feces_Shannon_5000_ASV_Y1
                  tab_spanner = c("**Specific richness**", "**Shannon diversity**"))}),
            group_header = spanner_names)

table_2_large <- 
  tbl_stack(tbls = list(
    tbl_merge(
      tbls = lapply(1:length(tbls_by_outcome_multi), function(i) tbls_by_outcome_multi[[i]][[1]]),  # correspond à ch_feces_SpecRich_5000_ASV_Y1
      tab_spanner = spanner_names),
    tbl_merge( 
      tbls = lapply(1:length(tbls_by_outcome_multi), function(i) tbls_by_outcome_multi[[i]][[2]]),  # correspond àch_feces_Shannon_5000_ASV_Y1
      tab_spanner = spanner_names)))

table_2 <- list(table_2_long = table_2_long, 
                table_2_large = table_2_large)
rm(table_2_long, table_2_large)

## Table 3: Associations betadiv and neurodev ----


## Table 4: Associations phyla and neurodev ----
table_4_long <-
  tbl_stack(tbls =
              lapply(tbls_by_outcome_multi, function(tbl_list) {
                tbl_merge(
                  list(tbl_list[[49]] %>% add_n(),   # correspond à ch_feces_rel_p1_Y1_10 (Firmicutes)
                       tbl_list[[50]] %>% add_n(),   # correspond à ch_feces_rel_p2_Y1_10 (Actinobacteria)
                       tbl_list[[51]] %>% add_n(),   # correspond à ch_feces_rel_p3_Y1_10 (Bacteroidetes)
                       tbl_list[[52]] %>% add_n()),  # correspond à ch_feces_rel_p4_Y1_10 (Proteobacteria)
                  tab_spanner = c("**Firmicutes**", 
                                  "**Actinobacteria**", 
                                  "**Bacteroidetes**", 
                                  "**Proteobacteria**"))}),
            group_header = spanner_names)

table_4_large <- 
  tbl_stack(tbls = list(
    tbl_merge(
      tbls = lapply(1:length(tbls_by_outcome_multi), function(i) tbls_by_outcome_multi[[i]][[49]]),  # correspond à ch_feces_rel_p1_Y1_10 (Firmicutes)
      tab_spanner = spanner_names),
    tbl_merge(
      tbls = lapply(1:length(tbls_by_outcome_multi), function(i) tbls_by_outcome_multi[[i]][[50]]),  # correspond à ch_feces_rel_p2_Y1_10 (Actinobacteria)
      tab_spanner = spanner_names),
    tbl_merge(
      tbls = lapply(1:length(tbls_by_outcome_multi), function(i) tbls_by_outcome_multi[[i]][[51]]),  # correspond à ch_feces_rel_p3_Y1_10 (Bacteroidetes)
      tab_spanner = spanner_names),
    tbl_merge(
      tbls = lapply(1:length(tbls_by_outcome_multi), function(i) tbls_by_outcome_multi[[i]][[52]]),  # correspond à ch_feces_rel_p4_Y1_10 (Proteobacteria)
      tab_spanner = spanner_names)))

table_4 <- list(table_4_long = table_4_long, 
                table_4_large = table_4_large)
rm(table_4_long, table_4_large)

# Figures ----
## Fig.1: Forestplot alpha div ----
figure_1 <- table_multi %>% 
  filter(`Gut microbiota parameters` %in% c("Shannon diversity",
                                             "Specific richness")) %>%
  mutate(Beta = as.numeric(Beta), 
         Outcome = 
           fct_relevel(Outcome,
                       "Total IQ score at 3 years", "Work memory IQ score at 3 years", 
                       "Visuospatiale IQ score at 3 years","Verbal comprehension IQ score at 3 years", 
                       "Plan and organization BRIEF-P score at 3 years", "Working memory BRIEF-P score at 3 years",
                       "Emotional control BRIEF-P score at 3 years","Shift BRIEF-P score at 3 years",
                       "Inhibition BRIEF-P score at 3 years", "Total SRS score at 3 years",
                       "Externalizing CBCL score at 2 years", "Internalizing CBCL score at 2 years"), 
         `Gut microbiota parameters` = 
           fct_relevel(`Gut microbiota parameters`, 
                       "Specific richness", "Shannon diversity")) %>%
  ggplot(aes(x = Outcome, 
             y = Beta, 
             min = lower_CI, 
             ymax = upper_CI, 
             color = p_value_shape)) +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_pointrange(
    position = position_dodge(width = 0.5), 
    size = 0.4) +
  labs(x = "Neurodevelopement", y = "", color = "p-value") +
  theme_lucid() +
  coord_flip()  +
  facet_wrap(~`Gut microbiota parameters`, scales = "free_x", ncol = 2) +
  theme(
    legend.position = "right",
    legend.box = "vertical", 
    legend.justification = "right") +
  scale_color_manual(values = c("p-value<0.05" = "red", "p-value≥0.05" = "black"))

figure_1
ggsave("4_output/fig.1 forest_plot_alphadiv.tiff", 
       figure_1, 
       device = "tiff",
       units = "cm",
       dpi = 300,
       height = 15, 
       width = 25)

## Fig.2: Boxplot beta div ----


## Fig.3: Forestplot phyla ----
figure_3 <- table_multi %>% 
  filter(`Gut microbiota parameters` %in% c("Firmicutes",
                                            "Actinobacteria", 
                                            "Bacteroidetes", 
                                            "Proteobacteria")) %>%
  mutate(Beta = as.numeric(Beta), 
         Outcome = 
           fct_relevel(Outcome,
                       "Total IQ score at 3 years", "Work memory IQ score at 3 years", 
                       "Visuospatiale IQ score at 3 years","Verbal comprehension IQ score at 3 years", 
                       "Plan and organization BRIEF-P score at 3 years", "Working memory BRIEF-P score at 3 years",
                       "Emotional control BRIEF-P score at 3 years","Shift BRIEF-P score at 3 years",
                       "Inhibition BRIEF-P score at 3 years", "Total SRS score at 3 years",
                       "Externalizing CBCL score at 2 years", "Internalizing CBCL score at 2 years"), 
         `Gut microbiota parameters` =
           fct_relevel(`Gut microbiota parameters`, 
                       "Firmicutes", "Actinobacteria", "Bacteroidetes", "Proteobacteria")) %>%
  ggplot(aes(x = Outcome, 
             y = Beta, 
             min = lower_CI, 
             ymax = upper_CI, 
             color = p_value_shape)) +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_pointrange(
    position = position_dodge(width = 0.5), 
    size = 0.4) +
  labs(x = "Neurodevelopement", y = "", color = "p-value") +
  theme_lucid() +
  coord_flip()  +
  facet_wrap(~`Gut microbiota parameters`, scales = "free_x", ncol = 4) +
  theme(
    legend.position = "right",
    legend.box = "vertical", 
    legend.justification = "right")  +
  scale_color_manual(values = c("p-value<0.05" = "red", "p-value≥0.05" = "black"))

figure_3
ggsave("4_output/fig.3 forest_plot_phyla.tiff", 
       figure_3, 
       device = "tiff",
       units = "cm",
       dpi = 300,
       height = 15, 
       width = 25)

## Fig.4: Mahatan plot genera ----
figure_4 <- table_multi  %>%
  filter(!`Gut microbiota parameters` %in% c("Firmicutes",
                                             "Actinobacteria",
                                             "Bacteroidetes",
                                             "Proteobacteria",
                                             "Shannon diversity",
                                             "Specific richness")) %>%
  ggplot(aes(x = -log10(`p-value`), y = `Gut microbiota parameters`)) +
  geom_point(aes(shape = sens_beta), size = 2) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = -log10(0.05/(33*7)), linetype = "dashed", color = "blue") +
  theme_lucid() +
  labs(x = "-log10(P-value)", 
       y = "Genera", 
       shape = "") +
  geom_text(aes(label = ifelse(`p-value` < 0.02, as.character(Outcome), "")), hjust = -0.05, vjust = -0.3, angle = 35, size = 3.5) +
  scale_shape_manual(values = c("Beta<0" = 15, "Beta≥0" = 17)) +# 15: carré plein, 17: triangle plein
  theme(
    legend.position = "right",
    legend.box = "vertical", 
    legend.justification = "right", 
    axis.text.y = element_text(face = "italic"))
figure_4
ggsave("4_output/fig.4 manhattan_plot_genera.tiff", 
       figure_4, 
       device = "tiff",
       units = "cm",
       dpi = 300,
       height = 25, 
       width = 40)

## Fig.5: Forestplot final genera ----
figure_5 <- table_multi %>% 
  filter(`p-value`<0.02) %>% 
  filter(!`Gut microbiota parameters` %in% c("Firmicutes",
                                             "Actinobacteria",
                                             "Bacteroidetes",
                                             "Proteobacteria",
                                             "Shannon diversity",
                                             "Specific richness")) %>%
  mutate(Beta = as.numeric(Beta), 
         Outcome = 
           fct_relevel(Outcome,
                       "Total IQ score at 3 years", "Work memory IQ score at 3 years", 
                       "Visuospatiale IQ score at 3 years","Verbal comprehension IQ score at 3 years", 
                       "Total SRS score at 3 years","Working memory BRIEF-P score at 3 years",
                       "Shift BRIEF-P score at 3 years","Plan and organization BRIEF-P score at 3 years",
                       "Inhibition BRIEF-P score at 3 years","Emotional control BRIEF-P score at 3 years",
                       "Externalizing CBCL score at 2 years", "Internalizing CBCL score at 2 years")) %>%
  ggplot(aes(x = Outcome, 
             y = Beta, 
             min = lower_CI, 
             ymax = upper_CI, 
             color = `Gut microbiota parameters`)) +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_pointrange(
    position = position_dodge(width = 0.5), 
    size = 0.4) +
  labs(x = "Neurodevelopement", y = "") +
  theme_lucid() +
  coord_flip()  +
  guides(color = guide_legend(title = "Genera"))+
  theme(
    legend.position = "right",
    legend.box = "vertical", 
    legend.justification = "right", 
    axis.text.y = element_text(face = "italic")) 

figure_5
ggsave("4_output/fig.5 forest_plot_genera.tiff", 
       figure_5, 
       device = "tiff",
       units = "cm",
       dpi = 300,
       height = 10, 
       width = 25)

# Additional tables ----
### Table S1: Distribution gut microbiota parameters ----
table_S1 <- descrip_num(data = bdd, 
                        vars = c("ch_feces_SpecRich_5000_ASV_Y1",
                                 "ch_feces_Shannon_5000_ASV_Y1",
                                 "ch_feces_rel_p1_Y1",         
                                 "ch_feces_rel_p2_Y1",
                                 "ch_feces_rel_p3_Y1",
                                 "ch_feces_rel_p4_Y1",
                                 genera_linear_complet))
write.xlsx(table_S1, file = "4_output/Table_S1.xlsx")

### Table S2: Distribution neurodevelopment ----
table_S2 <- descrip_num(data = bdd_final, vars = outcomes)
write.xlsx(table_S2, file = "4_output/Table_S2.xlsx")

### Table S3: Effects of the covariates on the outcomes ----
model_covars <- function(var_outcome, var_age, covars, data) {
  
  effectif <- data %>%                   # Création d'une colonne effectif
    filter(!is.na({{var_outcome}}))%>%
    select(all_of({{covars}}), 
           {{var_age}}) %>%
    tbl_summary(
      missing = "no", 
      type = list(ch_antibio_Y1_2cat ~ "categorical", 
                  home_total_y3 ~ "continuous"), 
      statistic = all_continuous() ~ "{N_nonmiss}") %>%
    bold_labels()
  
  tbl_model_univ <- data %>%                      # Création des colonnes modèle simple 
    select(
      {{var_outcome}}, 
      {{var_age}},
      all_of({{covars}})) %>%
    
    tbl_uvregression(
      method = lm ,
      y = {{var_outcome}},
      formula = "{y} ~ {x}", 
      hide_n = TRUE, 
      estimate_fun = scales::label_number(accuracy = .01, decimal.mark = "."),
      pvalue_fun = custom_pvalue_fun) %>%
    add_global_p() %>%
    bold_labels() %>%
    bold_p(t = 0.1)
  
  variables_explicatives <- c(covars, var_age)
  formule_regr <- as.formula(paste(var_outcome, "~", paste(variables_explicatives, collapse = " + ")))
  modele <- lm(formule_regr, data = data)
  tbl_model_multi <-  tbl_regression(
    modele, 
    estimate_fun = scales::label_number(accuracy = .01, decimal.mark = "."),
    pvalue_fun = custom_pvalue_fun) %>%
    add_global_p() %>%
    bold_labels() %>%
    bold_p(t = 0.1)
  
  tbl_model <- tbl_merge(
    tbls = list(effectif, 
                tbl_model_univ, 
                tbl_model_multi), 
    tab_spanner = c("", "**Univariate models**", "**Multivariate model**"))%>%
    modify_table_body(~.x %>% arrange(row_type == "glance_statistic"))
  
  return(tbl_model)
}

covariates <- c("po_w_kg_3cat",                                                 # redéfinition de covariates sinon il garde une variable age 
                "po_he_3cat",
                "mo_dipl_2cat",
                "mo_age",
                "mo_bmi_bepr_3cat",
                "ch_sex",
                "mo_par_2cat",
                "ch_bf_duration_till48w_4cat",
                "po_gd",
                "po_delmod",
                "ch_food_intro_Y1_3cat",
                "mo_pets",
                "ch_antibio_Y1_2cat",
                "home_total_y3",              
                "mo_hadtotscore_grt3_imp",
                "mo_tob_gr_anyt_yn_n2",
                "ch_tabacco_passive_up_to_Y1",
                "ch_care_main_12m_opt2_2c")

table_S3 <- tbl_merge(
  tbls = list(
    model_covars(var_outcome = "ch_cbclintscore_y2", var_age = "ch_age_CBCL_Y2", covars = covariates, data = bdd_final_imp_1),
    model_covars(var_outcome = "ch_cbclextscore_y2", var_age = "ch_age_CBCL_Y2", covars = covariates, data = bdd_final_imp_1),
    model_covars(var_outcome = "ch_SRStotal_y3", var_age = "ch_age_SRS_BRIEFP_Y3", covars = covariates, data = bdd_final_imp_1),
    model_covars(var_outcome = "ch_briefpinhibit_y3", var_age = "ch_age_SRS_BRIEFP_Y3", covars = covariates, data = bdd_final_imp_1),
    model_covars(var_outcome = "ch_briefpshift_y3", var_age = "ch_age_SRS_BRIEFP_Y3", covars = covariates, data = bdd_final_imp_1),
    model_covars(var_outcome = "ch_briefpemocontrol_y3", var_age = "ch_age_SRS_BRIEFP_Y3", covars = covariates, data = bdd_final_imp_1),
    model_covars(var_outcome = "ch_briefpworkmemo_y3", var_age = "ch_age_SRS_BRIEFP_Y3", covars = covariates, data = bdd_final_imp_1),
    model_covars(var_outcome = "ch_briefpplan_y3", var_age = "ch_age_SRS_BRIEFP_Y3", covars = covariates, data = bdd_final_imp_1),
    model_covars(var_outcome = "ch_verbal_comprehension_IQ_Y3", var_age = "ch_age_IQ_Y3", covars = covariates, data = bdd_final_imp_1),
    model_covars(var_outcome = "ch_visuospatiale_IQ_Y3", var_age = "ch_age_IQ_Y3", covars = covariates, data = bdd_final_imp_1),
    model_covars(var_outcome = "ch_work_memory_IQ_Y3", var_age = "ch_age_IQ_Y3", covars = covariates, data = bdd_final_imp_1), 
    model_covars(var_outcome = "ch_total_IQ_Y3", var_age = "ch_age_IQ_Y3", covars = covariates, data = bdd_final_imp_1)),
  tab_spanner = spanner_names)


## Table S4: Associations genera and neurodev ----
# Adjusted associations between the 46 most abundant genera in the child gut microbiota at one year and the neurodevelopment (n between X and X).
table_S4 <- tbl_stack(
  tbls = lapply(3:48, function(j) {  # tbl 3 à 48 : correspond aux analyses de genres
    tbl_merge(
      tbls = lapply(tbls_by_outcome_multi, function(i) i[[j]]),
      tab_spanner = spanner_names)}))

## Table S5: Sensitivity analysis - effect of the rarefaction threshold ----
# Sensitivity analysis - Adjusted associations between the gut microbiota α-diversity at different sequencing depths and the neurodevelopment. 
alpha_vec <- c("ch_feces_SpecRich_5000_ASV_Y1", "ch_feces_SpecRich_10000_ASV_Y1", 
               "ch_feces_Shannon_5000_ASV_Y1", "ch_feces_Shannon_10000_ASV_Y1")
table_S5 <- vector("list", length(outcomes))
names(table_S5) <- outcomes

# Boucle principale
for (outcome in outcomes) {
  # Sélection des covariables appropriées
  if (outcome %in% c("ch_cbclintscore_y2", "ch_cbclextscore_y2")) {
    covariates <- covariates_map$CBCL
  } else if (outcome %in% c("ch_SRStotal_y3", "ch_briefpinhibit_y3", "ch_briefpshift_y3", "ch_briefpemocontrol_y3", "ch_briefpworkmemo_y3", "ch_briefpplan_y3")) {
    covariates <- covariates_map$SRS_BRIEF
  } else if (outcome %in% c("ch_verbal_comprehension_IQ_Y3", "ch_visuospatiale_IQ_Y3", "ch_work_memory_IQ_Y3", "ch_total_IQ_Y3")) {
    covariates <- covariates_map$IQ
  }
  
  tbls_for_outcome_multi <- vector("list", length(alpha_vec))
  names(tbls_for_outcome_multi) <- alpha_vec
  
  for (exposure in alpha_vec) {                                               # running linear regression
    terms <- c(exposure, covariates)
    formula <- reformulate(terms, response = outcome)
    model <- lm(formula, data = bdd_final_imp_1_sensi_seuil)
    # model <- with(data = bdd_final_imp, 
    #               exp = lm(formula))
    
    tbl <-                                                                      
      tbl_regression(
        model, 
        include = exposure,
        estimate_fun = scales::label_number(accuracy = .01, decimal.mark = "."),
        pvalue_fun = custom_pvalue_fun,
        exponentiate = FALSE) %>%
      bold_p() %>%
      bold_labels() %>%
      add_global_p(include = exposure, singular.ok = TRUE, keep = TRUE) %>%
      add_n()%>%
      modify_table_body(~.x %>% arrange(row_type == "glance_statistic"))
    
    tbls_for_outcome_multi[[exposure]] <- tbl
  }
  table_S5[[outcome]] <- tbls_for_outcome_multi
}

table_S5_large <- tbl_stack(
  tbls = list(
    tbl_merge(tbls = lapply(tbls_by_outcome_multi, function(x) x[[1]]), tab_spanner = spanner_names),  # analyse principale seuil 5000
    tbl_merge(tbls = lapply(table_S5, function(x) x[[1]]), tab_spanner = spanner_names),               # analyse de sensibilité seuil 5000 pour n sur seuil 10000
    tbl_merge(tbls = lapply(table_S5, function(x) x[[2]]), tab_spanner = spanner_names),               # analyse de sensibilité seuil 10000 pour n sur seuil 10000
    tbl_merge(tbls = lapply(tbls_by_outcome_multi, function(x) x[[2]]), tab_spanner = spanner_names),  # analyse principale
    tbl_merge(tbls = lapply(table_S5, function(x) x[[3]]), tab_spanner = spanner_names),               # analyse de sensibilité seuil 5000 pour n sur seuil 10000
    tbl_merge(tbls = lapply(table_S5, function(x) x[[4]]), tab_spanner = spanner_names)))              # analyse de sensibilité seuil 10000 pour n sur seuil 10000

rm(terms, formula, model, exposure, outcome, alpha_vec, tbl)
covariates <- c("po_w_kg_3cat",                                                 # redéfinition de covariates sans variable age 
                "po_he_3cat",
                "mo_dipl_2cat",
                "mo_age",
                "mo_bmi_bepr_3cat",
                "ch_sex",
                "mo_par_2cat",
                "ch_bf_duration_till48w_4cat",
                "po_gd",
                "po_delmod",
                "ch_food_intro_Y1_3cat",
                "mo_pets",
                "ch_antibio_Y1_2cat",
                "home_total_y3",              
                "mo_hadtotscore_grt3_imp",
                "mo_tob_gr_anyt_yn_n2",
                "ch_tabacco_passive_up_to_Y1",
                "ch_care_main_12m_opt2_2c")

merge_tbls_function_rich <- function(index, tbls_by_outcome_multi, table_S5) {
  tbl_merge(
    tbls = list(
      tbls_by_outcome_multi[[index]][[1]] %>% add_n(), # Premier tbl_regression de tbls_by_outcome_multi
      table_S5[[index]][[1]],              # Premier tbl_regression de table_S5
      table_S5[[index]][[2]]               # Deuxième tbl_regression de table_S5
    ), 
    tab_spanner = c("**Threshold 5,000 (n=350)**","**Threshold 5,000 (n=339)**", "**Threshold 10,000 (n=339)**")
  )
}

merge_tbls_function_sha <- function(index, tbls_by_outcome_multi, table_S5) {
  tbl_merge(
    tbls = list(
      tbls_by_outcome_multi[[index]][[2]] %>% add_n(), # Premier tbl_regression de tbls_by_outcome_multi
      table_S5[[index]][[3]],              # Premier tbl_regression de table_S5
      table_S5[[index]][[4]]               # Deuxième tbl_regression de table_S5
    ), 
    tab_spanner = c("**Threshold 5,000 (n=350)**","**Threshold 5,000 (n=339)**", "**Threshold 10,000 (n=339)**")
  )
}

table_S5_long <- 
  tbl_stack(
    tbls = c(lapply(1:12, function(index) merge_tbls_function_rich(index, tbls_by_outcome_multi, table_S5)), 
             lapply(1:12, function(index) merge_tbls_function_sha(index, tbls_by_outcome_multi, table_S5))),
    group_header = c(spanner_names, spanner_names))
table_S5 <- list(table_S5_long = table_S5_long, 
                 table_S5_large = table_S5_large)
rm(merge_tbls_function_rich, merge_tbls_function_sha, table_S5_large, table_S5_long, tbls_for_outcome_multi)

## Table S6: Sensitivity analysis – effect of the HOME Y3 variable on CBCL Y2 -----
## Sensitivity analysis – Effects of the HOME covariate assessed at 3 years on the CBCL outcomes assessed at 2 years.
covariates_sensi_home <- 
  bdd_final_imp_1 %>%
  select(
    ch_age_CBCL_Y2,
    all_of(covariates)) %>%
  select(-home_total_y3) %>%
  colnames()
cbcl_vec <- c("ch_cbclintscore_y2", "ch_cbclextscore_y2")

# Initialisation de la liste pour stocker les résultats
prep_table_S6 <- list()

for (outcome in cbcl_vec) {
  # Initialisation de la sous-liste pour stocker les résultats pour chaque outcome
  table_S6_for_outcome <- list()
  
  for (explicative in explanatory) {
    # Construction de la formule
    formula <- as.formula(paste(outcome, "~", explicative, "+", paste(covariates_sensi_home, collapse = "+")))
    
    # Création du modèle de régression linéaire
    model <- lm(formula, data = bdd_final_imp_1)
    
    # Création du tbl_regression sans les covariates_sensi_home
    tbl <-                                                                      
      tbl_regression(
        model, 
        include = explicative,
        estimate_fun = scales::label_number(accuracy = .01, decimal.mark = "."),
        pvalue_fun = custom_pvalue_fun,
        exponentiate = FALSE) %>%
      bold_p() %>%
      bold_labels() %>%
      add_n()
    
    table_S6_for_outcome[[explicative]] <- tbl        # Stockage du tbl_regression dans la liste
  }
  
  prep_table_S6[[outcome]] <- table_S6_for_outcome    # Stockage des résultats pour chaque outcome
}
rm(cbcl_vec, outcome, explicative, formula, model, tbl, table_S6_for_outcome)


table_S6 <- list()

for (i in 1:52) {
  tbl <- tbl_merge(
    list(
      tbls_by_outcome_multi[[1]][[i]], # Premier tbl_regression de la première liste de "tbls_by_outcome_multi"
      prep_table_S6[[1]][[i]], # Premier tbl_regression de la première liste de "table_S6"
      tbls_by_outcome_multi[[2]][[i]], # Premier tbl_regression de la deuxième liste de "tbls_by_outcome_multi"
      prep_table_S6[[2]][[i]] # Premier tbl_regression de la deuxième liste de "table_S6"
    ),
    tab_spanner = c("**Internalizing CBCL score at 2 years, adjusted for HOME variable (main analysis)**", 
                    "**Internalizing CBCL score at 2 years, not adjusted for HOME variable**", 
                    "**Externalizing CBCL score at 2 years, adjusted for HOME variable (main analysis)**", 
                    "**Externalizing CBCL score at 2 years, not adjusted for HOME variable**")
  )
  table_S6[[i]] <- tbl   # Ajout de la table fusionnée à la liste
}
table_S6 <- tbl_stack(table_S6)

rm(covariates_sensi_home, i, tbl, prep_table_S6)


# Additional figures ----
## Figure S1: DAG ----
## Direct acyclic graph of the relation between one year child gut microbiota and neurodevelopment. 

## Figure S2: correlations gut microbiota parameters ----
## xx correlations between gut microbiota parameters (n between X and X). 
figure_S2 <- bdd %>% 
  select(all_of(genera_linear_complet)) 
colnames(figure_S2) <- genera_linear
figure_S2 <- figure_S2 %>%
  rename(
    "Escherichia and Shigella" = Escherichia_Shigella,
    "Clostridium XlVa" = Clostridium_XlVa, 
    "Lachnospiracea incertae sedis" = Lachnospiracea_incertae_sedis, 
    "Clostridium XVIII" = Clostridium_XVIII, 
    "Clostridium sensu stricto" = Clostridium_sensu_stricto, 
    "Ruminococcus 2" = Ruminococcus2, 
    "Clostridium IV" = Clostridium_IV, 
    "Erysipelotrichaceae incertae sedis" = Erysipelotrichaceae_incertae_sedis, 
    "Saccharibacteria genera incertae sedis" = Saccharibacteria_genera_incertae_sedis)

figure_S2 <- cor(figure_S2, 
                     use = "pairwise.complete.obs", 
                     method = "pearson")

plot.new()
tiff(filename = "4_output/Fig.S2 heatmap_cor_genera.tiff", units = "mm", width = 250, height = 250, res = 300)
corrplot(figure_S2, 
         method = 'color', 
         type = "lower", 
         tl.col = 'black', 
         tl.srt = 45, 
         # addCoef.col = "black",
         # number.cex = 0.5,
         # number.digits = 1,
         tl.cex = 0.5,
         col = rev(COL2(diverging = "RdYlBu")))
dev.off()

## Figure S3: correlations gut microbiota and covariates ----
## xx correlations between gut microbiota parameters and covariates (n between X and X).
figure_S3 <- bdd_final_imp_1 %>%
  select(all_of(outcomes),
         po_w_kg, 
         po_he, 
         mo_age, 
         mo_par,
         mo_bmi_bepr, 
         po_gd, 
         home_total_y3, 
         mo_hadtotscore_grt3_imp) %>%
  select_if(is.numeric)

figure_S3 <- round(cor(figure_S3,
                    use = "pairwise.complete.obs",
                    method = "pearson"), 1)

figure_S3 <- figure_S3 %>%
  as.data.frame() %>%
  select(!all_of(outcomes)) %>%
  t() %>%
  as.data.frame() %>%
  select(all_of(outcomes)) %>%
  as.matrix()
colnames(figure_S3) <- spanner_names

rownames(figure_S3) <- rownames(figure_S3) %>%
  str_replace_all(
    c("po_w_kg" = "Birth weight (kg)",
      "po_he" = "Birth length (cm)",
      "mo_age" = "Maternal age before pregnancy",
      "mo_par" = "Maternal parity",
      "mo_bmi_bepr" = "Maternal BMI before pregnancy",
      "po_gd" = "Gestational age (weeks)",
      "home_total_y3" = "HOME score at 3 years old",
      "mo_hadtotscore_grt3_imp" = "Maternal anxiety and depression during the third trim. of preg."))


tiff(filename = "4_output/Fig.S3 heatmap_cor_neuro_covar.tiff", units = "mm", width = 300, height = 300, res = 300)
corrplot(figure_S3,
         method = 'color',
         tl.col = 'black',
         tl.srt = 55,
         addCoef.col = "black",
         number.cex = 1,
         tl.cex = 1,
         number.digits = 1,
         col = rev(COL2(diverging = "RdYlBu")))
dev.off()


## Figure S4: correlations neurodevelopmental ----
## xx correlations between neurodevelopmental parameters (n between X and X).
figure_S4 <- bdd %>% 
  select(all_of(outcomes)) 
colnames(figure_S4) <- spanner_names
names(figure_S4) <- gsub(" score at 2 years", " Y2", names(figure_S4))
names(figure_S4) <- gsub(" score at 3 years", " Y3", names(figure_S4))

figure_S4 <- cor(figure_S4, 
                    use = "pairwise.complete.obs", 
                    method = "pearson")

plot.new()
tiff(filename = "4_output/Fig.S4 heatmap_cor_neuro.tiff", units = "mm", width = 250, height = 250, res = 300)
corrplot(figure_S4, 
         method = 'color', 
         type = "lower", 
         tl.col = 'black', 
         tl.srt = 45, 
         addCoef.col = "black",
         number.cex = 1,
         number.digits = 1,
         tl.cex = 1,
         col = rev(COL2(diverging = "RdYlBu")))
dev.off()

# Nettoyage final 
rm(covariates_CBCL, covariates_IQ, covariates_SRS_BRIEF, covariates_map)
# table_figure_list <- list(
#   table_1, 
#   table_2, 
#   # table_3,   # betadiv à faire
#   table_4, 
#   
#   figure_1, 
#   # table_2,   # betadiv à faire
#   figure_3, 
#   figure_4, 
#   figure_5, 
#   
#   table_S1, 
#   table_S2, 
#   table_S3, 
#   table_S4, 
#   table_S5,   
#   table_S6,  # sensi HOME CBCL à faire
#   
#   figure_S2, 
#   figure_S3, 
#   figure_S4)



# Strongest associations ----
## p<0.0035 ----
# 
# table_multi %>%
#   filter(is.na(Codage_rec)) %>%
#   filter(`p-value`<0.0035)%>%
#   View()
# 
# table_multi %>%
#   filter((Outcome == "Escherichia_Shigella" & Pollutants_Time_window_rec == "Benzophenone 3 M2") |
#            (Outcome == "Clostridium_XlVa" & Pollutants_Time_window_rec == "Benzophenone 3 t2")|
#            (Outcome == "Lachnospiracea_incertae_sedis" & Pollutants_Time_window_rec == "Butylparaben Y1") |
#            (Outcome == "Lachnospiracea_incertae_sedis" & Pollutants_Time_window_rec == "BPA Y1") | 
#            (Outcome == "Anaerostipes" & Pollutants_Time_window_rec == "Propylparaben t2") |
#            (Outcome == "Enterococcus" & Pollutants_Time_window_rec == "Ethylparaben t2") |
#            (Outcome == "Enterococcus" & Pollutants_Time_window_rec == "Butylparaben Y1") |
#            (Outcome == "Enterobacter" & Pollutants_Time_window_rec == "Ethylparaben Y1") |
#            (Outcome == "Collinsella" & Pollutants_Time_window_rec == "BPS M2") |
#            (Outcome == "Romboutsia" & Pollutants_Time_window_rec == "BPS M2") |
#            (Outcome == "Romboutsia" & Pollutants_Time_window_rec == "Benzophenone 3 t2") |
#            (Outcome == "Klebsiella" & Pollutants_Time_window_rec == "Ethylparaben Y1") |
#            (Outcome == "Coprococcus" & Pollutants_Time_window_rec == "BPA Y1") |
#            (Outcome == "Lactococcus" & Pollutants_Time_window_rec == "BPA M2") |
#            (Outcome == "Anaerotruncus" & Pollutants_Time_window_rec == "Butylparaben Y1")) %>%
#   View()
# 
# 
# results_signi <- tbl_merge(
#   tbls = list(
#     tbls_by_outcome_multi$Escherichia_Shigella$ch_OXBE_total_i_cor_M2_ln,
#     tbls_by_outcome_multi$Clostridium_XlVa$mo_OXBE_total_i_cor_t2_ln,
#     tbls_by_outcome_multi$Lachnospiracea_incertae_sedis$ch_BUPA_total_cat_Y1,
#     tbls_by_outcome_multi$Lachnospiracea_incertae_sedis$ch_BPA_total_i_cor_Y1_ln,
#     tbls_by_outcome_multi$Anaerostipes$mo_PRPA_total_i_cor_t2_ln,
#     tbls_by_outcome_multi$Enterococcus$mo_ETPA_total_i_cor_t2_ln,
#     tbls_by_outcome_multi$Enterococcus$ch_BUPA_total_cat_Y1,
#     tbls_by_outcome_multi$Enterobacter$ch_ETPA_total_i_cor_Y1_ln,
#     tbls_by_outcome_multi$Collinsella$ch_BPS_total_cat_M2_2,
#     tbls_by_outcome_multi$Romboutsia$ch_BPS_total_cat_M2_2,
#     tbls_by_outcome_multi$Romboutsia$mo_OXBE_total_i_cor_t2_ln,
#     tbls_by_outcome_multi$Klebsiella$ch_ETPA_total_i_cor_Y1_ln,
#     tbls_by_outcome_multi$Coprococcus$ch_BPA_total_i_cor_Y1_ln,
#     tbls_by_outcome_multi$Lactococcus$ch_BPA_total_i_cor_M2_ln,
#     tbls_by_outcome_multi$Anaerotruncus$ch_BUPA_total_cat_Y1), 
#   tab_spanner = c("**Escherichia_Shigella**", 
#                   "**Clostridium_XlVa**",
#                   "**Lachnospiracea_incertae_sedis**",
#                   "**Lachnospiracea_incertae_sedis**",
#                   "**Anaerostipes**",
#                   "**Enterococcus**",
#                   "**Enterococcus**",
#                   "**Enterobacter**",
#                   "**Collinsella**",
#                   "**Romboutsia**",
#                   "**Romboutsia**",
#                   "**Klebsiella**",
#                   "**Coprococcus**",
#                   "**Lactococcus**",
#                   "**Anaerotruncus**"))
# 
# write_xlsx(
#   x = taxa_table %>%
#     select(-"ch_feces_ASV_ID_Y1", -"ch_feces_domain_ASVbased_Y1", -"ch_feces_TAX_ASVbased_Y1") %>%
#     filter(ch_feces_genus_ASVbased_Y1 %in% c("Escherichia_Shigella", 
#                                              "Clostridium_XlVa",
#                                              "Lachnospiracea_incertae_sedis",
#                                              "Anaerostipes",
#                                              "Enterococcus",
#                                              "Enterobacter",
#                                              "Collinsella",
#                                              "Romboutsia",
#                                              "Klebsiella",
#                                              "Coprococcus",
#                                              "Lactococcus",
#                                              "Anaerotruncus")) %>% 
#     mutate_if(is.character, as.factor) %>%
#     arrange( ch_feces_phylum_ASVbased_Y1, 
#              ch_feces_class_ASVbased_Y1,
#              ch_feces_order_ASVbased_Y1,  
#              ch_feces_family_ASVbased_Y1 ) %>%
#     distinct(), 
#   path = "4_output/taxa manual/correspondnace taxa des principales asso.xlsx")
# 
# 
# ## p<0.001 ----
# table_multi %>%
#   filter(is.na(Codage_rec)) %>%
#   filter(`p-value`<0.001)%>%
#   View()
# 
# 
# table_multi %>%
#   filter((Outcome == "Escherichia_Shigella" & Pollutants_Time_window_rec == "Benzophenone 3 M2") |
#            (Outcome == "Clostridium_XlVa" & Pollutants_Time_window_rec == "Benzophenone 3 t2")|
#            (Outcome == "Lachnospiracea_incertae_sedis" & Pollutants_Time_window_rec == "Butylparaben Y1") |
#            (Outcome == "Lachnospiracea_incertae_sedis" & Pollutants_Time_window_rec == "BPA Y1") | 
#            (Outcome == "Anaerostipes" & Pollutants_Time_window_rec == "Propylparaben t2") |
#            (Outcome == "Enterococcus" & Pollutants_Time_window_rec == "Ethylparaben t2") |
#            (Outcome == "Enterococcus" & Pollutants_Time_window_rec == "Butylparaben Y1") |
#            (Outcome == "Enterobacter" & Pollutants_Time_window_rec == "Ethylparaben Y1") |
#            (Outcome == "Collinsella" & Pollutants_Time_window_rec == "BPS M2") |
#            (Outcome == "Romboutsia" & Pollutants_Time_window_rec == "BPS M2") |
#            (Outcome == "Romboutsia" & Pollutants_Time_window_rec == "Benzophenone 3 t2") |
#            (Outcome == "Klebsiella" & Pollutants_Time_window_rec == "Ethylparaben Y1") |
#            (Outcome == "Coprococcus" & Pollutants_Time_window_rec == "BPA Y1") |
#            (Outcome == "Lactococcus" & Pollutants_Time_window_rec == "BPA M2") |
#            (Outcome == "Anaerotruncus" & Pollutants_Time_window_rec == "Butylparaben Y1")) %>%
#   View()
# 
# 
# results_signi_0.001 <- tbl_merge(
#   tbls = list(
#     # tbls_by_outcome_multi$Escherichia_Shigella$ch_OXBE_total_i_cor_M2_ln,
#     # tbls_by_outcome_multi$Clostridium_XlVa$mo_OXBE_total_i_cor_t2_ln,
#     tbls_by_outcome_multi$Lachnospiracea_incertae_sedis$ch_BUPA_total_cat_Y1,
#     tbls_by_outcome_multi$Lachnospiracea_incertae_sedis$ch_BPA_total_i_cor_Y1_ln,
#     # tbls_by_outcome_multi$Anaerostipes$mo_PRPA_total_i_cor_t2_ln,
#     # tbls_by_outcome_multi$Enterococcus$mo_ETPA_total_i_cor_t2_ln,
#     tbls_by_outcome_multi$Enterococcus$ch_BUPA_total_cat_Y1,
#     tbls_by_outcome_multi$Enterobacter$ch_ETPA_total_i_cor_Y1_ln,
#     tbls_by_outcome_multi$Collinsella$ch_BPS_total_cat_M2_2,
#     # tbls_by_outcome_multi$Romboutsia$ch_BPS_total_cat_M2_2,
#     # tbls_by_outcome_multi$Romboutsia$mo_OXBE_total_i_cor_t2_ln,
#     tbls_by_outcome_multi$Klebsiella$ch_ETPA_total_i_cor_Y1_ln,
#     # tbls_by_outcome_multi$Coprococcus$ch_BPA_total_i_cor_Y1_ln,
#     tbls_by_outcome_multi$Lactococcus$ch_BPA_total_i_cor_M2_ln
#     # tbls_by_outcome_multi$Anaerotruncus$ch_BUPA_total_cat_Y1
#   ), 
#   tab_spanner = c(
#     # "**Escherichia_Shigella**", 
#     # "**Clostridium_XlVa**",
#     "**Lachnospiracea_incertae_sedis**",
#     "**Lachnospiracea_incertae_sedis**",
#     # "**Anaerostipes**",
#     # "**Enterococcus**",
#     "**Enterococcus**",
#     "**Enterobacter**",
#     "**Collinsella**",
#     # "**Romboutsia**",
#     # "**Romboutsia**",
#     "**Klebsiella**",
#     # "**Coprococcus**",
#     "**Lactococcus**"
#     # "**Anaerotruncus**"
#   ))
# 
# 
