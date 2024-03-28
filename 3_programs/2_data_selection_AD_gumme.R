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
    
    ch_WPPSI_verbal_comprehension_cor_Y3,
    ch_WPPSI_visuospatiale_cor_Y3,
    ch_WPPSI_work_memory_cor_Y3,
    ch_WPPSI_total_cor_Y3) %>%
  colnames()

spanner_names <- c("Internalizing CBCL score at 2 years",
                   "Externalizing CBCL score at 2 years",
                   "Total SRS score at 3 years",
                   "Inhibition BRIEF-P score at 3 years",
                   "Shift BRIEF-P score at 3 years",
                   "Emotional control BRIEF-P score at 3 years",
                   "Work memory BRIEF-P score at 3 years",
                   "Plan and organization BRIEF-P score at 3 years",
                   "Verbal comprehension WPPSI score at 3 years",
                   "Visuospatial WPPSI score at 3 years",
                   "Work memory WPPSI score at 3 years",
                   "Total WPPSI score at 3 years")


## Microbiote ----
microbiote_vec <- bdd %>% 
  select(
    ch_feces_SpecRich_5000_ASV_Y1_10,        
    ch_feces_Shannon_5000_ASV_Y1, 
    
    ch_feces_SpecRich_10000_ASV_Y1_10,      # pour analyses de sensibilité   
    ch_feces_Shannon_10000_ASV_Y1, 
    
    ch_feces_rel_p1_Y1_10, 
    ch_feces_rel_p2_Y1_10, 
    ch_feces_rel_p3_Y1_10, 
    ch_feces_rel_p4_Y1_10,
    
    all_of(genera_var_names)) %>%
  colnames()

## Covariables ----
covar_a_imputer <- bdd %>%
  select(
    ch_age_CBCL_Y2,
    ch_age_WPPSI_Y3,
    ch_age_SRS_BRIEFP_Y3, 
    ch_WPPSI_psy_Y3,
    
    
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
covar_vec_model_final <- bdd %>%
  select(
    # ne pas oublier d'ajouter covariable âge selon l'outcome de neuro (pour CBCL, SRS et BRIEF-P)
    # ne pas oublier d'ajouter covariable psy pour les outcomes WPPSI (en analyse de sensibilité ou principale à définir)
    
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
    home_total_y3,              
    mo_hadtotscore_grt3_imp,   
    mo_tob_gr_anyt_yn_n2,
    ch_tabacco_passive_up_to_Y1,
    ch_care_main_12m_opt2_2c) %>% 
  colnames()

# Export ----
save(list = ls(),
     file = "C:/Users/Aline/OneDrive - etu.univ-grenoble-alpes.fr/Documents/5. R projects/gut_microbiota_Y1_neurodev/1_intermediate_data/2_data_selection_AD_gumme.RData")
