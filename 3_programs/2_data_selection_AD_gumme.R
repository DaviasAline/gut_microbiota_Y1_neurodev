# 3_data_selection_AD_gumme
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

# Chargement des données 
load("1_intermediate_data/1_data_cleaning_AD_gumme.RData")
taxa_table_Y1 <- read_csv("0_source_data/taxa_table_ASVbased_Y1_AD_20220504_8.csv")  # ne pas inclure, table de correspondance taxonomique

# Création des vecteurs ----
## Neuro ----
neuro_vec <- bdd %>% 
  select(
    ch_cbclintscore_y2, 
    ch_cbclextscore_y2, 
    
    ch_socawar_y3, 
    ch_soccog_y3, 
    ch_soccom_y3, 
    ch_socmot_y3, 
    ch_RRB_y3, 
    ch_SRStotal_y3,
    
    ch_briefpinhibit_y3, 
    ch_briefpshift_y3,
    ch_briefpemocontrol_y3,
    ch_briefpworkmemo_y3, 
    ch_briefpplan_y3, 
    
    ch_verbal_comprehension_IQ_Y3,
    ch_visuospatiale_IQ_Y3,
    ch_work_memory_IQ_Y3,
    ch_total_IQ_Y3) %>%
  colnames()

## Microbiote ----
microbiote_vec <- bdd %>% 
  select(
    ch_feces_SpecRich_5000_ASV_Y1,        
    ch_feces_Shannon_5000_ASV_Y1, 
    
    ch_feces_SpecRich_10000_ASV_Y1,        
    ch_feces_Shannon_10000_ASV_Y1, 
    
    ch_feces_SpecRich_cmin_ASV_Y1, 
    ch_feces_Shannon_cmin_ASV_Y1,
    
    ch_feces_rel_p1_Y1, 
    ch_feces_rel_p2_Y1, 
    ch_feces_rel_p3_Y1, 
    ch_feces_rel_p4_Y1,
    
    all_of(genera_linear)) %>%
  colnames()

## Covariables ----
covar_a_imputer <- bdd %>%
  select(
    ch_age_CBCL_Y2,
    ch_age_IQ_Y3,
    ch_age_SRS_BRIEFP_Y3, 
    
    
    po_w_kg,      # utilisé précédement : po_w_kg_3cat
    po_w_kg_3cat,
    po_he,        # utilisé précédement : po_he_3cat
    po_he_3cat,
    ch_w_Y1,      # utilisé précédement : ch_w_Y1_3cat
    ch_w_Y1_3cat,
    ch_he_Y1,     # utilisé précédement : ch_he_Y1_3cat
    ch_he_Y1_3cat,
    ch_w_Y3,
    ch_w_Y3_3cat,
    ch_he_Y3,
    ch_he_Y3_3cat,
    ch_sex,
    
    mo_dipl_3cat,
    mo_dipl_2cat,
    mo_age,
    mo_age_3cat,
    mo_coupl_pr, 
    mo_bmi_bepr,  # utilisé précédement : mo_bmi_bepr_3cat
    mo_bmi_bepr_3cat, 
    mo_par,   # utilisé précédement : mo_par_2cat
    mo_par_2cat,
    ch_bf_duration_till48w,  # utilisé précédement : ch_bf_duration_till48w_4cat
    ch_bf_duration_till48w_4cat,
    po_gd, 
    po_delmod,
    ch_food_intro_Y1_3cat,
    mo_pets,
    ch_antibio_Y1_2cat,
    home_learning_y3,
    home_language_y3,
    home_academic_y3,
    home_variety_y3,
    home_total_y3, 
    mo_hadanxscore_grt3_imp, 
    mo_haddepscore_grt3_imp, 
    mo_hadtotscore_grt3_imp,
    mo_hadanxscore_y1_imp, 
    mo_haddepscore_y1_imp, 
    mo_hadtotscore_y1_imp, 
    mo_tob_gr_anyt_yn_n2,
    Mo_ETS_anyT_yn1_opt,
    ch_ETS_12m_opt36m,
    ch_tabacco_total_Y1,
    ch_tabacco_passive_up_to_Y1,
    mo_alcool_preg_opt_3cat,
    mo_alcool_preg_opt_4cat,
    ch_care_main_6m_opt2_2c, 
    ch_care_main_12m_opt2_2c, 
    ch_care_main_6m_12m_opt2_2c, 
    ch_care_main_6m_opt2_3c, 
    ch_care_main_12m_opt2_3c, 
    ch_care_main_6m_12m_opt2_3c,
    fa_ethnicity_2cat,  # exclue précédement 
    fa_dipl_3cat,
    fa_age,
    fa_age_3cat,
    mo_ethnicity_2cat,  # exclue précédement 
    ch_feces_age_w_Y1,
    ch_feces_RUN_Y1,
  ) %>%  
  colnames()

covar_age <- bdd %>%
  select(ch_age_CBCL_Y2,
         ch_age_IQ_Y3,
         ch_age_SRS_BRIEFP_Y3) %>%
  colnames()

covar_a_tester <- bdd %>%
  select(
    po_w_kg,      # utilisé précédement : po_w_kg_3cat
    po_w_kg_3cat,
    po_he,        # utilisé précédement : po_he_3cat
    po_he_3cat,
    ch_w_Y1,      # utilisé précédement : ch_w_Y1_3cat
    ch_w_Y1_3cat,
    ch_he_Y1,     # utilisé précédement : ch_he_Y1_3cat
    ch_he_Y1_3cat,
    ch_w_Y3,
    ch_w_Y3_3cat,
    ch_he_Y3,
    ch_he_Y3_3cat,
    ch_sex,
    
    mo_dipl_3cat,
    mo_dipl_2cat,
    mo_age,
    mo_bmi_bepr,  # utilisé précédement : mo_bmi_bepr_3cat
    mo_bmi_bepr_3cat, 
    mo_par_2cat,
    ch_bf_duration_till48w_4cat,
    po_gd, 
    po_delmod,
    ch_food_intro_Y1_3cat,
    mo_pets,
    ch_antibio_Y1_2cat,
    home_learning_y3,
    home_language_y3,
    home_academic_y3,
    home_variety_y3,
    home_total_y3, 
    mo_hadanxscore_grt3_imp, 
    mo_haddepscore_grt3_imp, 
    mo_hadtotscore_grt3_imp,
    mo_hadanxscore_y1_imp, 
    mo_haddepscore_y1_imp, 
    mo_hadtotscore_y1_imp, 
    mo_tob_gr_anyt_yn_n2,
    Mo_ETS_anyT_yn1_opt,
    ch_ETS_12m_opt36m,
    ch_tabacco_total_Y1,
    ch_tabacco_passive_up_to_Y1,
    mo_alcool_preg_opt_3cat,
    mo_alcool_preg_opt_4cat,
    ch_care_main_6m_opt2_2c, 
    ch_care_main_12m_opt2_2c, 
    ch_care_main_6m_12m_opt2_2c, 
    ch_care_main_6m_opt2_3c, 
    ch_care_main_12m_opt2_3c, 
    ch_care_main_6m_12m_opt2_3c,
    fa_dipl_3cat,
    fa_age,
    fa_age_3cat
    ) %>%  
  colnames()

## Vecteurs covariables
covar_vec_model_1 <- bdd %>%
  select(
    # ne pas oublier de changer covariable âge selon l'outcome de neuro 
    po_w_kg_3cat, 
    po_he_3cat, 
    mo_dipl_2cat,
    mo_age,
    mo_bmi_bepr_3cat,    
    ch_sex,
    mo_par_2cat, 
    ch_bf_duration_till48w_4cat,
    po_gd,
    po_delmod, 
    ch_food_intro_Y1_3cat,
    mo_pets,
    ch_antibio_Y1_2cat,
    home_total_y3,              # ne pas oublier de ne pas inclure cette variable dans le modèle CBCL 2 ans 
    mo_hadtotscore_grt3_imp,   
    mo_tob_gr_anyt_yn_n2,
    ch_tabacco_passive_up_to_Y1,
    
    ch_care_main_12m_opt2_2c) %>% 
  colnames()

# covar_vec_model_2 <- bdd %>%
#   select(
#     # ne pas oublier de changer selon l'outcome de neuro 
#     po_w_kg_3cat, 
#     po_he_3cat, 
#     mo_dipl_2cat,
#     mo_age,
#     mo_bmi_bepr_3cat,    
#     ch_sex,
#     mo_par_2cat, 
#     ch_bf_duration_till48w_4cat,
#     po_gd,
#     po_delmod, 
#     ch_food_intro_Y1_3cat,
#     mo_pets,
#     ch_antibio_Y1_2cat,
#     home_total_y3,            # ne pas oublier de ne pas inclure cette variable dans le modèle CBCL 2 ans 
#     mo_hadtotscore_grt3_imp,   
#     mo_tob_gr_anyt_yn_n2, 
#     ch_tabacco_passive_up_to_Y1,
#     
#     ch_care_main_6m_12m_opt2_2c) %>% 
#   colnames()

covar_vec_model_final <- covar_vec_model_1  # on choisit le modèle 1 (cf 5_stats)
rm(covar_vec_model_1)

# Export ----
save(list = ls(),
     file = "C:/Users/Aline/OneDrive - etu.univ-grenoble-alpes.fr/Documents/5. R projects/gut_microbiota_Y1_neurodev/1_intermediate_data/2_data_selection_AD_gumme.RData")
