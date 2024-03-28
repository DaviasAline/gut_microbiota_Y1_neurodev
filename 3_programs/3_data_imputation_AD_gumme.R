# 3_data_imputation_AD_gumme
# Aline Davias
# 25/08/2023

# Chargement des packages ----
library(readxl)
library(haven)
library(tidyverse)
library(lubridate)
library(gtsummary)
library(questionr)
library(mice)
library(expss)

# Chargement des données ----
load("1_intermediate_data/2_data_selection_AD_gumme.RData")
source("3_programs/5_functions_AD_gumme.R")

# Visualisation des données manquantes ----
bdd %>% 
  select(all_of(covar_a_imputer)) %>%
  tbl_summary(
    type = list(all_categorical() ~ "categorical", 
                home_learning_y3 ~ "continuous", 
                home_language_y3 ~ "continuous",    
                home_academic_y3 ~ "continuous",
                home_variety_y3 ~ "continuous",
                home_total_y3 ~ "continuous"), 
    digits = list(
      all_continuous() ~ 0,
      all_categorical() ~ c(0, 0)))%>% 
  bold_labels()

visu_na <- tbl_merge(
  tbls = list(
    tbl_1 = 
      bdd %>% 
      filter(!is.na(ch_cbclintscore_y2), !is.na(ch_feces_rel_p1_Y1)) %>%
      select(all_of(covar_a_imputer)) %>%
      select(-ch_age_WPPSI_Y3, -ch_age_SRS_BRIEFP_Y3, -ch_WPPSI_psy_Y3) %>%
      is.na() %>%
      as.data.frame() %>%
      tbl_summary(
        missing = "no",
        statistic = list(everything() ~ "{n} / {N} ({p}%)")) %>%
      bold_labels(), 
    tbl_2 = 
      bdd %>% 
      filter(!is.na(ch_cbclextscore_y2), !is.na(ch_feces_rel_p1_Y1)) %>%
      select(all_of(covar_a_imputer)) %>%
      select(-ch_age_WPPSI_Y3, -ch_age_SRS_BRIEFP_Y3, -ch_WPPSI_psy_Y3) %>%
      is.na() %>%
      as.data.frame() %>%
      tbl_summary(
        missing = "no",
        statistic = list(everything() ~ "{n} / {N} ({p}%)")) %>%
      bold_labels(), 
    tbl_3 = 
      bdd %>% 
      filter(!is.na(ch_SRStotal_y3), !is.na(ch_feces_rel_p1_Y1)) %>%
      select(all_of(covar_a_imputer)) %>%
      select(-ch_age_WPPSI_Y3, -ch_age_CBCL_Y2, -ch_WPPSI_psy_Y3) %>%
      is.na() %>%
      as.data.frame() %>%
      tbl_summary(
        missing = "no",
        statistic = list(everything() ~ "{n} / {N} ({p}%)")) %>%
      bold_labels(), 
    tbl_4 = 
      bdd %>% 
      filter(!is.na(ch_briefpinhibit_y3), !is.na(ch_feces_rel_p1_Y1)) %>%
      select(all_of(covar_a_imputer)) %>%
      select(-ch_age_WPPSI_Y3, -ch_age_CBCL_Y2, -ch_WPPSI_psy_Y3) %>%
      is.na() %>%
      as.data.frame() %>%
      tbl_summary(
        missing = "no",
        statistic = list(everything() ~ "{n} / {N} ({p}%)")) %>%
      bold_labels(), 
    tbl_5 = 
      bdd %>% 
      filter(!is.na(ch_briefpshift_y3), !is.na(ch_feces_rel_p1_Y1)) %>%
      select(all_of(covar_a_imputer)) %>%
      select(-ch_age_WPPSI_Y3, -ch_age_CBCL_Y2, -ch_WPPSI_psy_Y3) %>%
      is.na() %>%
      as.data.frame() %>%
      tbl_summary(
        missing = "no",
        statistic = list(everything() ~ "{n} / {N} ({p}%)")) %>%
      bold_labels(), 
    tbl_6 = 
      bdd %>% 
      filter(!is.na(ch_briefpemocontrol_y3), !is.na(ch_feces_rel_p1_Y1)) %>%
      select(all_of(covar_a_imputer)) %>%
      select(-ch_age_WPPSI_Y3, -ch_age_CBCL_Y2, -ch_WPPSI_psy_Y3) %>%
      is.na() %>%
      as.data.frame() %>%
      tbl_summary(
        missing = "no",
        statistic = list(everything() ~ "{n} / {N} ({p}%)")) %>%
      bold_labels(), 
    tbl_7 = 
      bdd %>% 
      filter(!is.na(ch_briefpworkmemo_y3), !is.na(ch_feces_rel_p1_Y1)) %>%
      select(all_of(covar_a_imputer)) %>%
      select(-ch_age_WPPSI_Y3, -ch_age_CBCL_Y2, -ch_WPPSI_psy_Y3) %>%
      is.na() %>%
      as.data.frame() %>%
      tbl_summary(
        missing = "no",
        statistic = list(everything() ~ "{n} / {N} ({p}%)")) %>%
      bold_labels(), 
    tbl_8 = 
      bdd %>% 
      filter(!is.na(ch_briefpplan_y3), !is.na(ch_feces_rel_p1_Y1)) %>%
      select(all_of(covar_a_imputer)) %>%
      select(-ch_age_WPPSI_Y3, -ch_age_CBCL_Y2, -ch_WPPSI_psy_Y3) %>%
      is.na() %>%
      as.data.frame() %>%
      tbl_summary(
        missing = "no",
        statistic = list(everything() ~ "{n} / {N} ({p}%)")) %>%
      bold_labels(), 
    tbl_9 = 
      bdd %>% 
      filter(!is.na(ch_WPPSI_verbal_comprehension_cor_Y3), !is.na(ch_feces_rel_p1_Y1)) %>%
      select(all_of(covar_a_imputer)) %>%
      select(-ch_age_SRS_BRIEFP_Y3, -ch_age_CBCL_Y2) %>%
      is.na() %>%
      as.data.frame() %>%
      tbl_summary(
        missing = "no",
        statistic = list(everything() ~ "{n} / {N} ({p}%)")) %>%
      bold_labels(), 
    tbl_10 = 
      bdd %>% 
      filter(!is.na(ch_WPPSI_visuospatiale_cor_Y3), !is.na(ch_feces_rel_p1_Y1)) %>%
      select(all_of(covar_a_imputer)) %>%
      select(-ch_age_SRS_BRIEFP_Y3, -ch_age_CBCL_Y2) %>%
      is.na() %>%
      as.data.frame() %>%
      tbl_summary(
        missing = "no",
        statistic = list(everything() ~ "{n} / {N} ({p}%)")) %>%
      bold_labels(), 
    tbl_11 = 
      bdd %>% 
      filter(!is.na(ch_WPPSI_work_memory_cor_Y3), !is.na(ch_feces_rel_p1_Y1)) %>%
      select(all_of(covar_a_imputer)) %>%
      select(-ch_age_SRS_BRIEFP_Y3, -ch_age_CBCL_Y2) %>%
      is.na() %>%
      as.data.frame() %>%
      tbl_summary(
        missing = "no",
        statistic = list(everything() ~ "{n} / {N} ({p}%)")) %>%
      bold_labels(), 
    tbl_12 = 
      bdd %>% 
      filter(!is.na(ch_WPPSI_total_cor_Y3), !is.na(ch_feces_rel_p1_Y1)) %>%
      select(all_of(covar_a_imputer)) %>%
      select(-ch_age_SRS_BRIEFP_Y3, -ch_age_CBCL_Y2) %>%
      is.na() %>%
      as.data.frame() %>%
      tbl_summary(
        missing = "no",
        statistic = list(everything() ~ "{n} / {N} ({p}%)")) %>%
      bold_labels()), 
  tab_spanner = c("**CBCL interne score Y2**", 
                  "**CBCL externe score Y2**", 
                  "**SRS total score Y3**", 
                  "**briefp inhibit y3**", 
                  "**briefp shift y3**", 
                  "**brief pemocontrol y3**", 
                  "**brief pworkmemo y3**", 
                  "**brief pplan y3**",
                  "**verbal comprehension WPPSI Y3**", 
                  "**visuospatial WPPSI Y3**", 
                  "**work memory WPPSI Y3**", 
                  "**total WPPSI Y3**"))
visu_na

# Imputation multiple des covariables ----
bdd_final <- bdd %>%
  select(ident, 
         all_of(neuro_vec), 
         all_of(microbiote_vec), 
         all_of(covar_a_imputer)) %>%
  filter(!is.na(ch_feces_rel_p1_Y1_10))

bdd_final <- bdd_final %>%
  mutate(
    fa_ethnicity_2cat = fct_recode(fa_ethnicity_2cat,
                                   "1" = "Caucassian",
                                   "2" = "Others"),
    mo_ethnicity_2cat = fct_recode(mo_ethnicity_2cat,
                                   "1" = "Caucassian",
                                   "2" = "Others"),
    mo_dipl_2cat = fct_recode(mo_dipl_2cat,
                              "1" = "≥5 years after graduation",
                              "2" = "<5 years after graduation"),
    po_delmod = fct_recode(po_delmod,
                           "1" = "Vaginal delivery",
                           "2" = "C-section"),
    mo_par_2cat = fct_recode(mo_par_2cat,
                             "1" = "None",
                             "2" = "1 child or more"),
    ch_sex = fct_recode(ch_sex,
                        "1" = "Male",
                        "2" = "Female"),
    mo_pets = fct_recode(mo_pets,
                         "1" = "No",
                         "2" = "One or more"),
    mo_tob_gr_anyt_yn_n2 = fct_recode(mo_tob_gr_anyt_yn_n2,
                                      "1" = "No",
                                      "2" = "Yes"),
    Mo_ETS_anyT_yn1_opt = fct_recode(Mo_ETS_anyT_yn1_opt,
                                     "1" = "No",
                                     "2" = "Yes"),

    ch_ETS_12m_opt36m = fct_recode(ch_ETS_12m_opt36m,
                                   "1" = "No",
                                   "2" = "Yes"),
    ch_tabacco_total_Y1 = fct_recode(ch_tabacco_total_Y1,
                                     "1" = "No",
                                     "2" = "Yes"),
    ch_tabacco_passive_up_to_Y1 = fct_recode(ch_tabacco_passive_up_to_Y1,
                                             "1" = "No",
                                             "2" = "Yes"),
    ch_antibio_Y1_2cat = fct_recode(ch_antibio_Y1_2cat,
                                    "1" = "No",
                                    "2" = "Yes"),
    ch_care_main_6m_opt2_2c = fct_recode(ch_care_main_6m_opt2_2c,
                                         "1" = "Collective day care",
                                         "2" = "Other"),
    ch_care_main_12m_opt2_2c = fct_recode(ch_care_main_12m_opt2_2c,
                                          "1" = "Collective day care",
                                          "2" = "Other"),
    ch_care_main_6m_12m_opt2_2c = fct_recode(ch_care_main_6m_12m_opt2_2c,
                                             "1" = "Collective day care",
                                             "2" = "Other"),
    ch_feces_RUN_Y1 = fct_recode(ch_feces_RUN_Y1,
                                 "1" = "R2",
                                 "2" = "R3"))

bdd_final_imp <- mice(bdd_final, 
                      seed = 11111,
                      maxit = 0)

method <- bdd_final_imp$method
method[c(microbiote_vec, neuro_vec)] <- ""
method

pred <- quickpred(bdd_final, 
                  exclude = neuro_vec, 
                  mincor = 0.2, 
                  minpuc = 0.4)


bdd_final_imp <- mice(bdd_final, m=5, meth = method, pred = pred, seed = 11111)

bdd_final_imp <- complete(bdd_final_imp, action = "long", include = TRUE)

bdd_final_imp <- bdd_final_imp %>%
  mutate(
    fa_ethnicity_2cat = fct_recode(fa_ethnicity_2cat,
                                   "Caucassian" = "1",
                                   "Others" = "2"),
    mo_ethnicity_2cat = fct_recode(mo_ethnicity_2cat,
                                   "Caucassian" = "1",
                                   "Others" = "2"),
    mo_dipl_2cat = fct_recode(mo_dipl_2cat,
                              "≥5 years after graduation" = "1",
                              "<5 years after graduation" = "2"),
    po_delmod = fct_recode(po_delmod,
                           "Vaginal delivery" = "1",
                           "C-section" = "2"),
    mo_par_2cat = fct_recode(mo_par_2cat,
                             "None" = "1",
                             "1 child or more" = "2"),
    ch_sex = fct_recode(ch_sex,
                        "Male" = "1",
                        "Female" = "2"),
    mo_pets = fct_recode(mo_pets,
                         "No" = "1",
                         "One or more" = "2"),
    mo_tob_gr_anyt_yn_n2 = fct_recode(mo_tob_gr_anyt_yn_n2,
                                      "No" = "1",
                                      "Yes" = "2"),
    Mo_ETS_anyT_yn1_opt = fct_recode(Mo_ETS_anyT_yn1_opt,
                                     "No" = "1",
                                     "Yes" = "2"),
    ch_ETS_12m_opt36m = fct_recode(ch_ETS_12m_opt36m,
                                   "No" = "1",
                                   "Yes" = "2"),
    ch_tabacco_total_Y1 = fct_recode(ch_tabacco_total_Y1,
                                     "No" = "1",
                                     "Yes" = "2"),
    ch_tabacco_passive_up_to_Y1 = fct_recode(ch_tabacco_passive_up_to_Y1,
                                             "No" = "1",
                                             "Yes" = "2"),
    ch_antibio_Y1_2cat = fct_recode(ch_antibio_Y1_2cat,
                                    "No" = "1",
                                    "Yes" = "2"),
    ch_care_main_6m_opt2_2c = fct_recode(ch_care_main_6m_opt2_2c,
                                         "Collective day care" = "1",
                                         "Other" = "2"),
    ch_care_main_12m_opt2_2c = fct_recode(ch_care_main_12m_opt2_2c,
                                          "Collective day care" = "1",
                                          "Other" = "2"),
    ch_care_main_6m_12m_opt2_2c = fct_recode(ch_care_main_6m_12m_opt2_2c,
                                             "Collective day care" = "1",
                                             "Other" = "2"),
    ch_feces_RUN_Y1 = fct_recode(ch_feces_RUN_Y1,
                                 "R2" = "1",
                                 "R3" = "2"))

bdd_final_imp_sensi_seuil <-                      # base de données réduite à n = 339 pour l'analyse de sensibilité du seuil de raréaction 
  bdd_final_imp %>%
  filter(!is.na(ch_feces_SpecRich_10000_ASV_Y1_10))
bdd_final_imp_sensi_seuil <- as.mids(bdd_final_imp_sensi_seuil)

bdd_final_imp <- as.mids(bdd_final_imp)

bdd_final <- bdd_final %>%
  mutate(
    fa_ethnicity_2cat = fct_recode(fa_ethnicity_2cat,
                                   "Caucassian" = "1",
                                   "Others" = "2"),
    mo_ethnicity_2cat = fct_recode(mo_ethnicity_2cat,
                                   "Caucassian" = "1",
                                   "Others" = "2"),
    mo_dipl_2cat = fct_recode(mo_dipl_2cat,
                              "≥5 years after graduation" = "1",
                              "<5 years after graduation" = "2"),
    po_delmod = fct_recode(po_delmod,
                           "Vaginal delivery" = "1",
                           "C-section" = "2"),
    mo_par_2cat = fct_recode(mo_par_2cat,
                             "None" = "1",
                             "1 child or more" = "2"),
    ch_sex = fct_recode(ch_sex,
                        "Male" = "1",
                        "Female" = "2"),
    mo_pets = fct_recode(mo_pets,
                         "No" = "1",
                         "One or more" = "2"),
    mo_tob_gr_anyt_yn_n2 = fct_recode(mo_tob_gr_anyt_yn_n2,
                                      "No" = "1",
                                      "Yes" = "2"),
    Mo_ETS_anyT_yn1_opt = fct_recode(Mo_ETS_anyT_yn1_opt,
                                     "No" = "1",
                                     "Yes" = "2"),
        ch_ETS_12m_opt36m = fct_recode(ch_ETS_12m_opt36m,
                                   "No" = "1",
                                   "Yes" = "2"),
    ch_tabacco_total_Y1 = fct_recode(ch_tabacco_total_Y1,
                                     "No" = "1",
                                     "Yes" = "2"),
    ch_tabacco_passive_up_to_Y1 = fct_recode(ch_tabacco_passive_up_to_Y1,
                                             "No" = "1",
                                             "Yes" = "2"),
    ch_antibio_Y1_2cat = fct_recode(ch_antibio_Y1_2cat,
                                    "No" = "1",
                                    "Yes" = "2"),
    ch_care_main_6m_opt2_2c = fct_recode(ch_care_main_6m_opt2_2c,
                                         "Collective day care" = "1",
                                         "Other" = "2"),
    ch_care_main_12m_opt2_2c = fct_recode(ch_care_main_12m_opt2_2c,
                                          "Collective day care" = "1",
                                          "Other" = "2"),
    ch_care_main_6m_12m_opt2_2c = fct_recode(ch_care_main_6m_12m_opt2_2c,
                                             "Collective day care" = "1",
                                             "Other" = "2"),
    ch_feces_RUN_Y1 = fct_recode(ch_feces_RUN_Y1,
                                 "R2" = "1",
                                 "R3" = "2"))

# Pour accéder aux différents jeux de données imputés :
bdd_final_imp_1 <- complete(bdd_final_imp, 1) # 1er jeu de données imputé
bdd_final_imp_2 <- complete(bdd_final_imp, 2) # 2eme jeu de données imputé
bdd_final_imp_3 <- complete(bdd_final_imp, 3) # 3eme jeu de données imputé
bdd_final_imp_4 <- complete(bdd_final_imp, 4) # 4eme jeu de données imputé
bdd_final_imp_5 <- complete(bdd_final_imp, 5) # 5eme jeu de données imputé

bdd_final_imp_1_sensi_seuil <- complete(bdd_final_imp_sensi_seuil, 1) # 1er jeu de données imputé
bdd_final_imp_2_sensi_seuil <- complete(bdd_final_imp_sensi_seuil, 2) # 2eme jeu de données imputé
bdd_final_imp_3_sensi_seuil <- complete(bdd_final_imp_sensi_seuil, 3) # 3eme jeu de données imputé
bdd_final_imp_4_sensi_seuil <- complete(bdd_final_imp_sensi_seuil, 4) # 4eme jeu de données imputé
bdd_final_imp_5_sensi_seuil <- complete(bdd_final_imp_sensi_seuil, 5) # 5eme jeu de données imputé


rm(pred, 
   #bdd_final_imp_1, 
   bdd_final_imp_2, bdd_final_imp_3, bdd_final_imp_4, bdd_final_imp_5, 
   method, 
   #bdd_final_imp_1_sensi_seuil, 
   bdd_final_imp_2_sensi_seuil, bdd_final_imp_3_sensi_seuil, bdd_final_imp_4_sensi_seuil, bdd_final_imp_5_sensi_seuil)  


# Vérification de l'imputation ----
# /!\ à revoir avec Alicia et Matthieu si il y a la possibilité de mettre include = neuro_vec sans laisser plein de NA ----
visu_na <- tbl_merge(
  tbls = list(
    tbl_1 = 
      bdd_final_imp_1 %>% 
      select(all_of(covar_a_imputer), ch_cbclintscore_y2) %>%
      select(-ch_age_WPPSI_Y3, -ch_age_SRS_BRIEFP_Y3, -ch_WPPSI_psy_Y3) %>%
      is.na() %>%
      as.data.frame() %>%
      tbl_summary(
        missing = "no",
        statistic = list(everything() ~ "{n} / {N} ({p}%)")) %>%
      bold_labels(), 
    tbl_2 = 
      bdd_final_imp_1 %>% 
      select(all_of(covar_a_imputer), ch_cbclextscore_y2) %>%
      select(-ch_age_WPPSI_Y3, -ch_age_SRS_BRIEFP_Y3, -ch_WPPSI_psy_Y3) %>%
      is.na() %>%
      as.data.frame() %>%
      tbl_summary(
        missing = "no",
        statistic = list(everything() ~ "{n} / {N} ({p}%)")) %>%
      bold_labels(), 
    tbl_3 = 
      bdd_final_imp_1 %>% 
      select(all_of(covar_a_imputer), ch_SRStotal_y3) %>%
      select(-ch_age_WPPSI_Y3, -ch_age_CBCL_Y2, -ch_WPPSI_psy_Y3) %>%
      is.na() %>%
      as.data.frame() %>%
      tbl_summary(
        missing = "no",
        statistic = list(everything() ~ "{n} / {N} ({p}%)")) %>%
      bold_labels(), 
    tbl_4 = 
      bdd_final_imp_1 %>% 
      select(all_of(covar_a_imputer), ch_briefpinhibit_y3) %>%
      select(-ch_age_WPPSI_Y3, -ch_age_CBCL_Y2, -ch_WPPSI_psy_Y3) %>%
      is.na() %>%
      as.data.frame() %>%
      tbl_summary(
        missing = "no",
        statistic = list(everything() ~ "{n} / {N} ({p}%)")) %>%
      bold_labels(), 
    tbl_5 = 
      bdd_final_imp_1 %>% 
      select(all_of(covar_a_imputer), ch_briefpshift_y3) %>%
      select(-ch_age_WPPSI_Y3, -ch_age_CBCL_Y2, -ch_WPPSI_psy_Y3) %>%
      is.na() %>%
      as.data.frame() %>%
      tbl_summary(
        missing = "no",
        statistic = list(everything() ~ "{n} / {N} ({p}%)")) %>%
      bold_labels(), 
    tbl_6 = 
      bdd_final_imp_1 %>% 
      select(all_of(covar_a_imputer), ch_briefpemocontrol_y3) %>%
      select(-ch_age_WPPSI_Y3, -ch_age_CBCL_Y2, -ch_WPPSI_psy_Y3) %>%
      is.na() %>%
      as.data.frame() %>%
      tbl_summary(
        missing = "no",
        statistic = list(everything() ~ "{n} / {N} ({p}%)")) %>%
      bold_labels(), 
    tbl_7 = 
      bdd_final_imp_1 %>% 
      select(all_of(covar_a_imputer), ch_briefpworkmemo_y3) %>%
      select(-ch_age_WPPSI_Y3, -ch_age_CBCL_Y2, -ch_WPPSI_psy_Y3) %>%
      is.na() %>%
      as.data.frame() %>%
      tbl_summary(
        missing = "no",
        statistic = list(everything() ~ "{n} / {N} ({p}%)")) %>%
      bold_labels(), 
    tbl_8 = 
      bdd_final_imp_1 %>% 
      select(all_of(covar_a_imputer), ch_briefpplan_y3) %>%
      select(-ch_age_WPPSI_Y3, -ch_age_CBCL_Y2, -ch_WPPSI_psy_Y3) %>%
      is.na() %>%
      as.data.frame() %>%
      tbl_summary(
        missing = "no",
        statistic = list(everything() ~ "{n} / {N} ({p}%)")) %>%
      bold_labels(), 
    tbl_9 = 
      bdd_final_imp_1 %>% 
      select(all_of(covar_a_imputer), ch_WPPSI_verbal_comprehension_cor_Y3) %>%
      select(-ch_age_SRS_BRIEFP_Y3, -ch_age_CBCL_Y2) %>%
      is.na() %>%
      as.data.frame() %>%
      tbl_summary(
        missing = "no",
        statistic = list(everything() ~ "{n} / {N} ({p}%)")) %>%
      bold_labels(), 
    tbl_10 = 
      bdd_final_imp_1 %>% 
      select(all_of(covar_a_imputer), ch_WPPSI_visuospatiale_cor_Y3) %>%
      select(-ch_age_SRS_BRIEFP_Y3, -ch_age_CBCL_Y2) %>%
      is.na() %>%
      as.data.frame() %>%
      tbl_summary(
        missing = "no",
        statistic = list(everything() ~ "{n} / {N} ({p}%)")) %>%
      bold_labels(), 
    tbl_11 = 
      bdd_final_imp_1 %>% 
      select(all_of(covar_a_imputer), ch_WPPSI_work_memory_cor_Y3) %>%
      select(-ch_age_SRS_BRIEFP_Y3, -ch_age_CBCL_Y2) %>%
      is.na() %>%
      as.data.frame() %>%
      tbl_summary(
        missing = "no",
        statistic = list(everything() ~ "{n} / {N} ({p}%)")) %>%
      bold_labels(), 
    tbl_12 = 
      bdd_final_imp_1 %>% 
      select(all_of(covar_a_imputer), ch_WPPSI_total_cor_Y3) %>%
      select(-ch_age_SRS_BRIEFP_Y3, -ch_age_CBCL_Y2) %>%
      is.na() %>%
      as.data.frame() %>%
      tbl_summary(
        missing = "no",
        statistic = list(everything() ~ "{n} / {N} ({p}%)")) %>%
      bold_labels()), 
  tab_spanner = c("**CBCL interne score Y2**", 
                  "**CBCL externe score Y2**", 
                  "**SRS total score Y3**", 
                  "**briefp inhibit y3**", 
                  "**briefp shift y3**", 
                  "**brief pemocontrol y3**", 
                  "**brief pworkmemo y3**", 
                  "**brief pplan y3**",
                  "**verbal comprehension WPPSI Y3**", 
                  "**visuospatial WPPSI Y3**", 
                  "**work memory WPPSI Y3**", 
                  "**total WPPSI Y3**"))
visu_na
rm(visu_na)

# Export ----
save(list = ls(),
     file = "1_intermediate_data/3_data_imputation_AD_gumme.RData")
