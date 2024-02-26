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
taxa_table <- read_csv("0_source_data/taxa_table_ASVbased_Y1_AD_20220504_8.csv")

# Vectors ----
covariates <- covar_vec_model_final
explanatory <- bdd %>% 
  select(all_of(microbiote_vec)) %>% 
  select(-ch_feces_SpecRich_cmin_ASV_Y1, 
         -ch_feces_Shannon_cmin_ASV_Y1,
         -ch_feces_SpecRich_10000_ASV_Y1, 
         -ch_feces_Shannon_10000_ASV_Y1) %>%
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


### Heatmap of correlation ----
cormat_genera <- bdd %>% 
  select(all_of(genera_linear_complet)) 
colnames(cormat_genera) <- genera_linear
cormat_genera <- cormat_genera %>%
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

cormat_genera <- cor(cormat_genera, 
                           use = "pairwise.complete.obs", 
                           method = "pearson")

corrplot(cormat_genera, 
         method = 'color', 
         type = "lower", 
         tl.col = 'black', 
         tl.srt = 45, 
         # addCoef.col = "black",
         # number.cex = 0.5,
         # number.digits = 1,
         tl.cex = 0.6,
         col = rev(COL2(diverging = "RdYlBu")))


# Calcul du nombre de tests ----



# Running linear regressions ----
corres <- 
  taxa_table %>% 
  select(Phyla_corres = ch_feces_phylum_ASVbased_Y1, 
         Class_corres = ch_feces_class_ASVbased_Y1,
         Order_corres = ch_feces_order_ASVbased_Y1, 
         Family_corres = ch_feces_family_ASVbased_Y1,
         Exposure = ch_feces_genus_ASVbased_Y1) %>%
  filter(Exposure %in% genera_linear) %>%
  distinct(Exposure, .keep_all = TRUE)

covariates_CBCL <- c("ch_age_CBCL_Y2", covariates)
covariates_IQ <- c("ch_age_IQ_Y3", covariates)
covariates_SRS_BRIEF <- c("ch_age_SRS_BRIEFP_Y3", covariates)

## Version multivariée ----
####### Test ######
covariates_map <- list(
  CBCL = covariates_CBCL,
  SRS_BRIEF = covariates_SRS_BRIEF,
  IQ = covariates_IQ
)

# Création d'une liste pour stocker les tableaux par outcome
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
  
  for (exposure in explanatory) {
    terms <- c(exposure, covariates)
    formula <- reformulate(terms, response = outcome)
    model <- lm(formula, data = bdd_final_imp_1)
    # model <- with(data = bdd_final_imp, 
    #               exp = lm(formula))
    
    tbl <-                                                                      # running linear regression
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


### Tableau pour l'article ----
stacked_tbls_by_outcome_multi <- vector("list", length(outcomes))            # Création liste pour stocker les tableaux empilés par outcome
names(stacked_tbls_by_outcome_multi) <- outcomes

for (outcome in names(tbls_by_outcome_multi)) {                                   # Récupérer les tableaux de régression pour cet outcome
  tbls_for_this_outcome_multi <- tbls_by_outcome_multi[[outcome]]
  stacked_tbl_multi <- do.call(tbl_stack, list(tbls = tbls_for_this_outcome_multi)) # Empiler les tableaux en un seul tableau
  stacked_tbls_by_outcome_multi[[outcome]] <- stacked_tbl_multi                     # Ajouter le tableau empilé à la liste des tableaux empilés
}
labels_outcomes <- bdd_final_imp_1 %>%
  select(all_of(outcomes)) %>%
  sapply(function(x) attr(x, "label"))
results_tbl_multi <- tbl_merge(tbls = stacked_tbls_by_outcome_multi,                # Fusionner les tableaux empilés en un seul tableau
                             tab_spanner = labels_outcomes)


### Tableau pour générer des figures ----
table_multi <- tibble()                                                           # Initialisation d'un tibble vide pour stocker les résultats finaux
for (i in seq_along(tbls_by_outcome_multi)) {                                     # Nom de l'outcome pour cette itération
  outcome_name <- names(tbls_by_outcome_multi)[i]
  
  for (j in seq_along(tbls_by_outcome_multi[[i]])) {                              # Itération sur chaque tbl_regression dans la liste courante
    exposure_name <- names(tbls_by_outcome_multi[[i]])[j]                         # Nom de la variable d'exposition pour cette itération
    tbl_data <- tbls_by_outcome_multi[[i]][[j]] %>%                               # Extraction des données du tableau tbl_regression
      as_tibble() %>%
      mutate(Outcome = outcome_name, Exposure = exposure_name)
    
    table_multi <- bind_rows(table_multi, tbl_data)                                 # Ajout des données extraites au tibble final
  }
}
rm(terms, formula, model, 
   tbl, tbl_data, tbls_for_this_outcome_multi, 
   stacked_tbl_multi, tbls_for_outcome_multi,
   stacked_tbls_by_outcome_multi,
   exposure, exposure_name, i, j, outcome, outcome_name)

test <-                                                                    # Ajout variable de la correspondance en phyla
  left_join(table_multi, corres, by = "Exposure") %>%
  select(Phyla_corres, Class_corres, Order_corres, Family_corres, everything())

test <- test %>%
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
                 "Firmicutes" = "ch_feces_rel_p1_Y1",
                 "Actinobacteria" = "ch_feces_rel_p2_Y1",
                 "Bacteroidetes" = "ch_feces_rel_p3_Y1",
                 "Proteobacteria" = "ch_feces_rel_p4_Y1",
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
    `q-value` = `p-value`/(29*31), 
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
descrip_covar

## Table 2: Associations alphadiv and neurodev ----


## Table 3: Associations betadiv and neurodev ----


## Table 4: Associations phyla and neurodev ----


# Figures ----
## Fig.1: Forestplot alpha div ----
forest_plot_alphadiv <- test %>% 
  filter(`Gut microbiota parameters` %in% c("Shannon diversity",
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
    legend.justification = "right") 

forest_plot_alphadiv
ggsave("4_output/fig.1 forest_plot_alphadiv.tiff", 
       forest_plot_alphadiv, 
       device = "tiff",
       units = "cm",
       dpi = 300,
       height = 15, 
       width = 25)

## Fig.3: Forestplot phyla ----
forest_plot_phyla <- test %>% 
  filter(`Gut microbiota parameters` %in% c("Firmicutes",
                                            "Actinobacteria", 
                                            "Bacteroidetes", 
                                            "Proteobacteria")) %>%
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
    legend.justification = "right") 

forest_plot_phyla
ggsave("4_output/fig.3 forest_plot_phyla.tiff", 
       forest_plot_phyla, 
       device = "tiff",
       units = "cm",
       dpi = 300,
       height = 15, 
       width = 25)

## Fig.4: Mahatan plot genera ----
mahatan_plot <- test  %>%
  filter(!`Gut microbiota parameters` %in% c("Firmicutes",
                                             "Actinobacteria",
                                             "Bacteroidetes",
                                             "Proteobacteria",
                                             "Shannon diversity",
                                             "Specific richness")) %>%
  ggplot(aes(x = -log10(`p-value`), y = `Gut microbiota parameters`)) +
  geom_point(aes(shape = sens_beta), size = 2) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = -log10(0.05/(29*31)), linetype = "dashed", color = "blue") +
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
mahatan_plot
ggsave("4_output/fig.4 manhattan_plot_genera.tiff", 
       mahatan_plot, 
       device = "tiff",
       units = "cm",
       dpi = 300,
       height = 20, 
       width = 40)

## Fig.5: Forestplot final genera ----
forest_plot_genera <- test %>% 
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

forest_plot_genera
ggsave("4_output/fig.5 forest_plot_genera.tiff", 
       forest_plot_genera, 
       device = "tiff",
       units = "cm",
       dpi = 300,
       height = 15, 
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

### Table S2: distribution neurodevelopment ----
table_S2 <- descrip_num(data = bdd, vars = outcomes)

# Additional figures ----


# Strongest associations ----
## p<0.0035 ----

table_multi %>%
  filter(is.na(Codage_rec)) %>%
  filter(`p-value`<0.0035)%>%
  View()

table_multi %>%
  filter((Outcome == "Escherichia_Shigella" & Pollutants_Time_window_rec == "Benzophenone 3 M2") |
           (Outcome == "Clostridium_XlVa" & Pollutants_Time_window_rec == "Benzophenone 3 t2")|
           (Outcome == "Lachnospiracea_incertae_sedis" & Pollutants_Time_window_rec == "Butylparaben Y1") |
           (Outcome == "Lachnospiracea_incertae_sedis" & Pollutants_Time_window_rec == "BPA Y1") | 
           (Outcome == "Anaerostipes" & Pollutants_Time_window_rec == "Propylparaben t2") |
           (Outcome == "Enterococcus" & Pollutants_Time_window_rec == "Ethylparaben t2") |
           (Outcome == "Enterococcus" & Pollutants_Time_window_rec == "Butylparaben Y1") |
           (Outcome == "Enterobacter" & Pollutants_Time_window_rec == "Ethylparaben Y1") |
           (Outcome == "Collinsella" & Pollutants_Time_window_rec == "BPS M2") |
           (Outcome == "Romboutsia" & Pollutants_Time_window_rec == "BPS M2") |
           (Outcome == "Romboutsia" & Pollutants_Time_window_rec == "Benzophenone 3 t2") |
           (Outcome == "Klebsiella" & Pollutants_Time_window_rec == "Ethylparaben Y1") |
           (Outcome == "Coprococcus" & Pollutants_Time_window_rec == "BPA Y1") |
           (Outcome == "Lactococcus" & Pollutants_Time_window_rec == "BPA M2") |
           (Outcome == "Anaerotruncus" & Pollutants_Time_window_rec == "Butylparaben Y1")) %>%
  View()


results_signi <- tbl_merge(
  tbls = list(
    tbls_by_outcome_multi$Escherichia_Shigella$ch_OXBE_total_i_cor_M2_ln,
    tbls_by_outcome_multi$Clostridium_XlVa$mo_OXBE_total_i_cor_t2_ln,
    tbls_by_outcome_multi$Lachnospiracea_incertae_sedis$ch_BUPA_total_cat_Y1,
    tbls_by_outcome_multi$Lachnospiracea_incertae_sedis$ch_BPA_total_i_cor_Y1_ln,
    tbls_by_outcome_multi$Anaerostipes$mo_PRPA_total_i_cor_t2_ln,
    tbls_by_outcome_multi$Enterococcus$mo_ETPA_total_i_cor_t2_ln,
    tbls_by_outcome_multi$Enterococcus$ch_BUPA_total_cat_Y1,
    tbls_by_outcome_multi$Enterobacter$ch_ETPA_total_i_cor_Y1_ln,
    tbls_by_outcome_multi$Collinsella$ch_BPS_total_cat_M2_2,
    tbls_by_outcome_multi$Romboutsia$ch_BPS_total_cat_M2_2,
    tbls_by_outcome_multi$Romboutsia$mo_OXBE_total_i_cor_t2_ln,
    tbls_by_outcome_multi$Klebsiella$ch_ETPA_total_i_cor_Y1_ln,
    tbls_by_outcome_multi$Coprococcus$ch_BPA_total_i_cor_Y1_ln,
    tbls_by_outcome_multi$Lactococcus$ch_BPA_total_i_cor_M2_ln,
    tbls_by_outcome_multi$Anaerotruncus$ch_BUPA_total_cat_Y1), 
  tab_spanner = c("**Escherichia_Shigella**", 
                  "**Clostridium_XlVa**",
                  "**Lachnospiracea_incertae_sedis**",
                  "**Lachnospiracea_incertae_sedis**",
                  "**Anaerostipes**",
                  "**Enterococcus**",
                  "**Enterococcus**",
                  "**Enterobacter**",
                  "**Collinsella**",
                  "**Romboutsia**",
                  "**Romboutsia**",
                  "**Klebsiella**",
                  "**Coprococcus**",
                  "**Lactococcus**",
                  "**Anaerotruncus**"))

write_xlsx(
  x = taxa_table %>%
    select(-"ch_feces_ASV_ID_Y1", -"ch_feces_domain_ASVbased_Y1", -"ch_feces_TAX_ASVbased_Y1") %>%
    filter(ch_feces_genus_ASVbased_Y1 %in% c("Escherichia_Shigella", 
                                             "Clostridium_XlVa",
                                             "Lachnospiracea_incertae_sedis",
                                             "Anaerostipes",
                                             "Enterococcus",
                                             "Enterobacter",
                                             "Collinsella",
                                             "Romboutsia",
                                             "Klebsiella",
                                             "Coprococcus",
                                             "Lactococcus",
                                             "Anaerotruncus")) %>% 
    mutate_if(is.character, as.factor) %>%
    arrange( ch_feces_phylum_ASVbased_Y1, 
             ch_feces_class_ASVbased_Y1,
             ch_feces_order_ASVbased_Y1,  
             ch_feces_family_ASVbased_Y1 ) %>%
    distinct(), 
  path = "4_output/taxa manual/correspondnace taxa des principales asso.xlsx")


## p<0.001 ----
table_multi %>%
  filter(is.na(Codage_rec)) %>%
  filter(`p-value`<0.001)%>%
  View()


table_multi %>%
  filter((Outcome == "Escherichia_Shigella" & Pollutants_Time_window_rec == "Benzophenone 3 M2") |
           (Outcome == "Clostridium_XlVa" & Pollutants_Time_window_rec == "Benzophenone 3 t2")|
           (Outcome == "Lachnospiracea_incertae_sedis" & Pollutants_Time_window_rec == "Butylparaben Y1") |
           (Outcome == "Lachnospiracea_incertae_sedis" & Pollutants_Time_window_rec == "BPA Y1") | 
           (Outcome == "Anaerostipes" & Pollutants_Time_window_rec == "Propylparaben t2") |
           (Outcome == "Enterococcus" & Pollutants_Time_window_rec == "Ethylparaben t2") |
           (Outcome == "Enterococcus" & Pollutants_Time_window_rec == "Butylparaben Y1") |
           (Outcome == "Enterobacter" & Pollutants_Time_window_rec == "Ethylparaben Y1") |
           (Outcome == "Collinsella" & Pollutants_Time_window_rec == "BPS M2") |
           (Outcome == "Romboutsia" & Pollutants_Time_window_rec == "BPS M2") |
           (Outcome == "Romboutsia" & Pollutants_Time_window_rec == "Benzophenone 3 t2") |
           (Outcome == "Klebsiella" & Pollutants_Time_window_rec == "Ethylparaben Y1") |
           (Outcome == "Coprococcus" & Pollutants_Time_window_rec == "BPA Y1") |
           (Outcome == "Lactococcus" & Pollutants_Time_window_rec == "BPA M2") |
           (Outcome == "Anaerotruncus" & Pollutants_Time_window_rec == "Butylparaben Y1")) %>%
  View()


results_signi_0.001 <- tbl_merge(
  tbls = list(
    # tbls_by_outcome_multi$Escherichia_Shigella$ch_OXBE_total_i_cor_M2_ln,
    # tbls_by_outcome_multi$Clostridium_XlVa$mo_OXBE_total_i_cor_t2_ln,
    tbls_by_outcome_multi$Lachnospiracea_incertae_sedis$ch_BUPA_total_cat_Y1,
    tbls_by_outcome_multi$Lachnospiracea_incertae_sedis$ch_BPA_total_i_cor_Y1_ln,
    # tbls_by_outcome_multi$Anaerostipes$mo_PRPA_total_i_cor_t2_ln,
    # tbls_by_outcome_multi$Enterococcus$mo_ETPA_total_i_cor_t2_ln,
    tbls_by_outcome_multi$Enterococcus$ch_BUPA_total_cat_Y1,
    tbls_by_outcome_multi$Enterobacter$ch_ETPA_total_i_cor_Y1_ln,
    tbls_by_outcome_multi$Collinsella$ch_BPS_total_cat_M2_2,
    # tbls_by_outcome_multi$Romboutsia$ch_BPS_total_cat_M2_2,
    # tbls_by_outcome_multi$Romboutsia$mo_OXBE_total_i_cor_t2_ln,
    tbls_by_outcome_multi$Klebsiella$ch_ETPA_total_i_cor_Y1_ln,
    # tbls_by_outcome_multi$Coprococcus$ch_BPA_total_i_cor_Y1_ln,
    tbls_by_outcome_multi$Lactococcus$ch_BPA_total_i_cor_M2_ln
    # tbls_by_outcome_multi$Anaerotruncus$ch_BUPA_total_cat_Y1
  ), 
  tab_spanner = c(
    # "**Escherichia_Shigella**", 
    # "**Clostridium_XlVa**",
    "**Lachnospiracea_incertae_sedis**",
    "**Lachnospiracea_incertae_sedis**",
    # "**Anaerostipes**",
    # "**Enterococcus**",
    "**Enterococcus**",
    "**Enterobacter**",
    "**Collinsella**",
    # "**Romboutsia**",
    # "**Romboutsia**",
    "**Klebsiella**",
    # "**Coprococcus**",
    "**Lactococcus**"
    # "**Anaerotruncus**"
  ))


