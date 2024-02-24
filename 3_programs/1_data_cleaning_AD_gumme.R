# 1_data_cleaning_AD_gumme
# Aline Davias
# 25/08/2023

# Chargement des packages ----
library(haven)
library(tidyverse)
library(lubridate)
library(gtsummary)
library(questionr)
library(mice)
library(expss)

# Chargement des données ----
load("C:/Users/Aline/OneDrive - etu.univ-grenoble-alpes.fr/Documents/5. R projects/gut_microbiota_Y1_neurodev/1_intermediate_data/0_source_data_reading_AD_gumme.RData")
taxa_table_Y1 <- read_csv("C:/Users/Aline/OneDrive - etu.univ-grenoble-alpes.fr/Documents/5. R projects/gut_microbiota_Y1_neurodev/0_source_data/taxa_table_ASVbased_Y1_AD_20220504_8.csv")  # ne pas inclure, table de correspondance taxonomique

# Nettoyage des données ----
## Covariates ----
bdd <- bdd %>%
  rename(
    fa_ethnicity_2cat = fa_ethnicity,
    mo_ethnicity_2cat = mo_ethnicity, 
    mo_ethnicity = mt1saa1_q02,
    ch_sibling = mt1saa1_q04,
    fa_dipl = ft2sac1_q03,
    mo_pets = mt3eaf1_q01,       # Variables animaux (questionnaire MT3EAF1_V1)
    mo_dogs = mt3eaf1_q01p1,
    mo_cats = mt3eaf1_q01p5,
    mo_birds = mt3eaf1_q01p9,
    mo_rodents = mt3eaf1_q01p13,
    mo_other_pets = mt3eaf1_q01p17, 
    ch_verbal_comprehension_IQ_Y3 = cy3hax1_q15, 
    ch_visuospatiale_IQ_Y3 = cy3hax1_q16, 
    ch_work_memory_IQ_Y3 = cy3hax1_q17, 
    ch_total_IQ_Y3 = cy3hax1_q18, 
    ch_date_IQ_Y3 = cy3hax1_q19)


bdd[, c("fa_ethnicity_2cat",
        "fa_dipl",
        
        "mo_ethnicity_2cat",
        "mo_ethnicity",
        "mo_dipl",
        "mo_coupl_pr", 
        "mo_tob_gr_anyt_yn_n2",
        "Mo_ETS_anyT_yn1_opt",
        "mo_alcool_preg",
        "mo_alcool_preg_opt", 
        
        "ch_sex",
        "mo_par",
        "ch_ETS_12m_opt36m",
        "po_delmod",
        "ch_food_intro_Y1",
        "mo_pets",
        "ch_feces_RUN_Y1", 
        "ch_antibio_Y1")] <- lapply(bdd[, c("fa_ethnicity_2cat",
                                       "fa_dipl",
                                       
                                       "mo_ethnicity_2cat",
                                       "mo_ethnicity",
                                       "mo_dipl",
                                       "mo_coupl_pr", 
                                       "mo_tob_gr_anyt_yn_n2",
                                       "Mo_ETS_anyT_yn1_opt",
                                       "mo_alcool_preg",
                                       "mo_alcool_preg_opt", 
                                       
                                       "ch_sex",
                                       "mo_par",
                                       "ch_ETS_12m_opt36m",
                                       "po_delmod",
                                       "ch_food_intro_Y1",
                                       "mo_pets",
                                       "ch_feces_RUN_Y1", 
                                       "ch_antibio_Y1")], as.character)
bdd[, c("fa_ethnicity_2cat",
        "fa_dipl",
        
        "mo_ethnicity_2cat",
        "mo_ethnicity",
        "mo_dipl",
        "mo_coupl_pr", 
        "mo_tob_gr_anyt_yn_n2",
        "Mo_ETS_anyT_yn1_opt",
        "mo_alcool_preg",
        "mo_alcool_preg_opt", 
        
        "ch_sex",
        "mo_par",
        "ch_ETS_12m_opt36m",
        "po_delmod",
        "ch_food_intro_Y1",
        "mo_pets",
        "ch_feces_RUN_Y1", 
        "ch_antibio_Y1")] <- lapply(bdd[, c("fa_ethnicity_2cat",
                                       "fa_dipl",
                                       
                                       "mo_ethnicity_2cat",
                                       "mo_ethnicity",
                                       "mo_dipl",
                                       "mo_coupl_pr", 
                                       "mo_tob_gr_anyt_yn_n2",
                                       "Mo_ETS_anyT_yn1_opt",
                                       "mo_alcool_preg",
                                       "mo_alcool_preg_opt", 
                                       
                                       "ch_sex",
                                       "mo_par",
                                       "ch_ETS_12m_opt36m",
                                       "po_delmod",
                                       "ch_food_intro_Y1",
                                       "mo_pets",
                                       "ch_feces_RUN_Y1", 
                                       "ch_antibio_Y1")], as.factor)

bdd <- bdd %>%
  mutate(
    mo_tob_gr_anyt_yn_n2 = fct_recode(mo_tob_gr_anyt_yn_n2,            # tabagisme actif de la mère pdt la grossesse
                                      "No" = "0",
                                      "Yes" = "1"), 
    ch_ETS_12m_opt36m = fct_recode(ch_ETS_12m_opt36m,                  # tabagisme passif de l'enfant 
                                   "No" = "0",
                                   "Yes" = "1"), 
    Mo_ETS_anyT_yn1_opt = fct_recode(Mo_ETS_anyT_yn1_opt,              # tabagisme passif de la mère pdt la grossesse
                                     "No" = "0",
                                     "Yes" = "1"), 
    ch_sex = fct_recode(ch_sex,                                        # Labellisation de ch_sex
                        "Male" = "1", 
                        "Female" = "2"), 
    po_delmod = fct_recode(po_delmod,                                  # Labellisation de po_delmod
                           "Vaginal delivery" = "1", 
                           "C-section" = "2"), 
    mo_par_2cat = fct_recode(mo_par,                                        # Labellisation de mo_par
                             "None" = "0", 
                             "1 child or more" = "1", 
                             "1 child or more" = "2"), 
    mo_ethnicity = fct_recode(mo_ethnicity,                            # Labellisation de mo_ethnicity
                              "Afrique" = "1",
                              "Amériques" = "2",
                              "Asie du Sud-Est" = "3",
                              "Europe" = "4",
                              "Méditéranée Orientale" = "5",
                              "Pacifique Occidental" = "6",            
                              "Autre" = "7", 
                              "Ne sait pas/ne souhaite pas répondre" = "99"), 
    mo_dipl = fct_recode(mo_dipl,                                      # Labellisation des variables éducation mère mo_dipl
                         "BEP/CAP/Highschool" = "2", 
                         "1-2 years after graduation" = "3", 
                         "3-4 years after graduation" = "4", 
                         "≥5 years after graduation" = "5"), 
    fa_dipl = fct_recode(fa_dipl,                                      # Labellisation des variables éducation père fa_dipl              
                         "BEP/CAP/Highschool" = "1", 
                         "BEP/CAP/Highschool" = "2", 
                         "BEP/CAP/Highschool" = "3", 
                         "BEP/CAP/Highschool" = "4", 
                         "BEP/CAP/Highschool" = "5", 
                         "1-2 years after graduation" = "6",
                         "3-4 years after graduation" = "7",
                         "≥5 years after graduation" = "8"), 
    mo_pets = fct_recode(mo_pets,                                      # Labellisation des variables animaux 
                         "No" = "0", 
                         "One or more" = "1"), 
    mo_alcool_preg = fct_recode(
      mo_alcool_preg,
      "No drinking" = "1",
      "<1drink/month when she didn’t know she was pregnant and did not drink when she knew she was pregnant (or did not answer)" = "2",
      ">1drink/month when she didn’t know she was pregnant and did not drink when she knew she was pregnant" = "3",
      "<1drink/month when she knew she was pregnant" = "4",
      ">1drink/month when she knew she was pregnant" = "5"), 
    mo_alcool_preg_opt = fct_recode(
      mo_alcool_preg_opt,
      "No drinking" = "1",
      "<1drink/month when she didn’t know she was pregnant and did not drink when she knew she was pregnant (or did not answer)" = "2",
      ">1drink/month when she didn’t know she was pregnant and did not drink when she knew she was pregnant" = "3",
      "<1drink/month when she knew she was pregnant" = "4",
      ">1drink/month when she knew she was pregnant" = "5"), 
    fa_dipl = fct_relevel(fa_dipl, 
                          "≥5 years after graduation", 
                          "3-4 years after graduation", 
                          "1-2 years after graduation",
                          "BEP/CAP/Highschool"))
# warning message mo_ethnicity : Unknown levels in `f`: 6  : normal car absence d'observation pour cette modalité 
# warning message fa_dipl : Unknown levels in `f`: 1, 2  : normal car absence d'observation pour ces modalités 


## Gut microbiota ----
## Imputation des 0 par la valeur 1/5000 + log transformation des variables genres 
bdd_genera_imp <- bdd %>%
  select(ident, 
         contains("ch_feces_rel_g")) %>%
  column_to_rownames("ident")
var_label(bdd_genera_imp) <-  str_replace(
  var_label(bdd_genera_imp),                                    # set correct variable names                   
  "One year child feces relative abundance of genus ", "")
colnames(bdd_genera_imp) <- var_label(bdd_genera_imp)

genera_linear <- bdd_genera_imp %>%                             # create a vector 
  na.omit() %>%
  select_if(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.3) %>%
  colnames()

bdd_genera_imp <- bdd_genera_imp %>% 
  select(all_of(genera_linear)) %>%                     # selection des genres ayant une abondance >30%
  mutate_all(., ~ ifelse(. == 0, 1/5000, .)) %>%        # remplacement des valeurs 0 par 1/5000
  mutate_all(~ log(.)) %>%                              # transformation logarithmique
  rename_with(~gsub("genus ", "", .), everything()) %>% # changement des noms de colonnes pour qu'ils n'aient pas d'espace
  rownames_to_column(var = "ident")




# Création variable âge du père ----
bdd <- bdd %>% 
  mutate(
    fa_age = (difftime(po_datedel, fa_bd, units = "days"))/365, 
    fa_age = as.numeric(fa_age))



# Création des variables poids et taille de l'enfant à 3 ans ----
bdd <- bdd %>%
  mutate(
    ch_w_Y3 = case_when(
    !is.na(cy3cak1_q04) & !is.na(cy3cak1_q05) ~ (cy3cak1_q04 + cy3cak1_q05)/2, 
    is.na(cy3cak1_q04) & !is.na(cy3cak1_q05) ~ cy3cak1_q05, 
    !is.na(cy3cak1_q04) & is.na(cy3cak1_q05) ~ cy3cak1_q04, 
    .default = NA), 
    ch_he_Y3 = case_when(
      ident == "27309" ~ NA,    # valeurs improbables, attendre réponse Sarah ou Anne 
      ident == "17668" ~ NA, # valeurs improbables, attendre réponse Sarah ou Anne 
      ident == "17958" ~ NA, # valeurs improbables, attendre réponse Sarah ou Anne 
      !is.na(cy3cak1_q06) & !is.na(cy3cak1_q07) ~ (cy3cak1_q06 + cy3cak1_q07)/2, 
      is.na(cy3cak1_q06) & !is.na(cy3cak1_q07) ~ cy3cak1_q07, 
      !is.na(cy3cak1_q06) & is.na(cy3cak1_q07) ~ cy3cak1_q06, 
      .default = NA))


# Création des variables âges de l'enfant aux tests de neurodeveloppement  ----
bdd <- bdd %>%
  mutate(
    ch_age_SRS_BRIEFP_Y3 = difftime(cy3hat1_date_creation, po_datedel, units = "weeks"), 
    ch_age_SRS_BRIEFP_Y3 = as.numeric(ch_age_SRS_BRIEFP_Y3),
    ch_age_CBCL_Y2 = difftime(cy2ebc1_date_creation, po_datedel, units = "weeks"),
    ch_age_CBCL_Y2 = as.numeric(ch_age_CBCL_Y2),
    ch_age_IQ_Y3 = difftime(cy3hax1_date_creation, po_datedel, units = "weeks"),     # Le test du QI a été réalisé en même temsp que l'examen clinique des 3 ans 
    ch_age_IQ_Y3 = as.numeric(ch_age_IQ_Y3),
    ch_age_Vineland_Y2 = difftime(cy2hap2_date_creation, po_datedel, units = "weeks"), 
    ch_age_Vineland_Y2 = as.numeric(ch_age_Vineland_Y2),
    ch_age_Vineland_Y2 = ifelse(ident == "27201", NA, ch_age_Vineland_Y2),
    ch_age_Vineland_Y2_bis = difftime(cy2hap3_date_creation, po_datedel, units = "weeks"), 
    ch_age_Vineland_Y2_bis = as.numeric(ch_age_Vineland_Y2_bis))


# Création d'une variable exposition périnatale au tabac ----
bdd <- bdd %>%
  mutate(
    ch_tabacco_total_Y1 = case_when(mo_tob_gr_anyt_yn_n2 == "Yes" | Mo_ETS_anyT_yn1_opt == "Yes" | ch_ETS_12m_opt36m == "Yes" ~ "Yes",  
                                    mo_tob_gr_anyt_yn_n2 == "No" & Mo_ETS_anyT_yn1_opt == "No" & ch_ETS_12m_opt36m == "No" ~ "No",   
                                    
                                    is.na(mo_tob_gr_anyt_yn_n2) & is.na(Mo_ETS_anyT_yn1_opt) & is.na(ch_ETS_12m_opt36m) ~ NA,
                                    
                                    is.na(mo_tob_gr_anyt_yn_n2) & Mo_ETS_anyT_yn1_opt == "No" & ch_ETS_12m_opt36m == "No" ~ NA,
                                    mo_tob_gr_anyt_yn_n2 == "No" & is.na(Mo_ETS_anyT_yn1_opt) & ch_ETS_12m_opt36m == "No" ~ NA,
                                    mo_tob_gr_anyt_yn_n2 == "No" & Mo_ETS_anyT_yn1_opt == "No" & is.na(ch_ETS_12m_opt36m) ~ NA,
                                    
                                    is.na(mo_tob_gr_anyt_yn_n2) & is.na(Mo_ETS_anyT_yn1_opt) & ch_ETS_12m_opt36m == "No" ~ NA,
                                    is.na(mo_tob_gr_anyt_yn_n2) & Mo_ETS_anyT_yn1_opt == "No" & is.na(ch_ETS_12m_opt36m) ~ NA,
                                    mo_tob_gr_anyt_yn_n2 == "No" & is.na(Mo_ETS_anyT_yn1_opt) & is.na(ch_ETS_12m_opt36m) ~ NA), 
    
    ch_tabacco_passive_up_to_Y1 = case_when(Mo_ETS_anyT_yn1_opt == "Yes" | ch_ETS_12m_opt36m == "Yes" ~ "Yes",  
                                           Mo_ETS_anyT_yn1_opt == "No" & ch_ETS_12m_opt36m == "No" ~ "No",   
                                           
                                           is.na(Mo_ETS_anyT_yn1_opt) & is.na(ch_ETS_12m_opt36m) ~ NA,
                                           
                                           is.na(Mo_ETS_anyT_yn1_opt) & ch_ETS_12m_opt36m == "No" ~ NA,
                                           Mo_ETS_anyT_yn1_opt == "No" & is.na(ch_ETS_12m_opt36m) ~ NA), 
    ch_tabacco_total_Y1 = as.factor(ch_tabacco_total_Y1), 
    ch_tabacco_passive_up_to_Y1 = as.factor(ch_tabacco_passive_up_to_Y1)) 


# Création variable mode de garde up to one year ----
bdd <- bdd %>% 
  mutate(
    ch_care_main_6m_12m_opt2_2c = 
      case_when(ch_care_main_6m_opt2_2c == "Collective day care" & ch_care_main_12m_opt2_2c == "Collective day care" ~ "Collective day care", 
                ch_care_main_6m_opt2_2c == "Collective day care" | ch_care_main_12m_opt2_2c == "Collective day care" ~ "Collective day care", 
                ch_care_main_6m_opt2_2c == "Other" & ch_care_main_12m_opt2_2c == "Other" ~ "Other", 
                ch_care_main_6m_opt2_2c == "Other" | ch_care_main_12m_opt2_2c == "Other" ~ "Other", 
                .default = NA), 
    ch_care_main_6m_12m_opt2_3c = 
      case_when(ch_care_main_6m_opt2_3c == "Collective day care" & ch_care_main_12m_opt2_3c == "Collective day care" ~ "Collective day care", 
                ch_care_main_6m_opt2_3c == "Collective day care" | ch_care_main_12m_opt2_3c == "Collective day care" ~ "Collective day care", 
                ch_care_main_6m_opt2_3c == "Nursery assistant (out child home)" & ch_care_main_12m_opt2_3c == "Nursery assistant (out child home)" ~ "Nursery assistant (out child home)", 
                ch_care_main_6m_opt2_3c == "Nursery assistant (out child home)" | ch_care_main_12m_opt2_3c == "Nursery assistant (out child home)" ~ "Nursery assistant (out child home)", 
                ch_care_main_6m_opt2_3c == "Other" & ch_care_main_12m_opt2_3c == "Other" ~ "Other", 
                ch_care_main_6m_opt2_3c == "Other" | ch_care_main_12m_opt2_3c == "Other" ~ "Other", 
                .default = NA), 
    ch_care_main_6m_12m_opt2_2c = as.factor(ch_care_main_6m_12m_opt2_2c), 
    ch_care_main_6m_12m_opt2_3c = as.factor(ch_care_main_6m_12m_opt2_3c))

# Choix du codage des covariables ----
bdd <- bdd %>%                                   
  mutate(
    
    
    mo_ethnicity_2cat = fct_recode(mo_ethnicity_2cat,                           # ethinicité de la mère
                                   "Caucassian" = "1",
                                   "Others" = "2"), 
    fa_ethnicity_2cat = fct_recode(fa_ethnicity_2cat,                           # ethinicité du père
                                   "Caucassian" = "1", 
                                   "Others" = "2"), 
    
    mo_coupl_pr = fct_recode(mo_coupl_pr,                                       # statut matrimonial
                             "No" = "0",
                             "Yes" = "1"), 
    
    mo_age_4cat = cut(mo_age,                                                   # âge de la mère    
                      include.lowest = TRUE,
                      right = FALSE,
                      dig.lab = 4,
                      breaks = c(20, 27, 33, 36, 46)), 
    mo_age_4cat = fct_recode(mo_age_4cat,                                      
                             "20-26 years" = "[20,27)",
                             "27-32 years" = "[27,33)",
                             "33-35 years" = "[33,36)",
                             ">35 years" = "[36,46]"), 
    mo_age_3cat = cut(mo_age,
                      include.lowest = TRUE,
                      right = FALSE,
                      dig.lab = 1,
                      breaks = c(20, 30.5598907470703, 34.064338684082, 45.4702262878418)), 
    mo_age_3cat = fct_recode(mo_age_3cat,
                             "20-30 years" = "[20,31)",
                             "31-33 years" = "[31,34)",
                             "34-44 years" = "[34,45]"), 
    fa_age_3cat = cut(fa_age,                                                   # âge du père
                      include.lowest = TRUE,
                      right = FALSE,
                      dig.lab = 2,
                      breaks = c(23, 32.4977168949772, 36.7013698630137, 73)),
    fa_age_3cat = fct_recode(fa_age_3cat,
                             "23-31 years" = "[23,32)",
                             "32-36 years" = "[32,37)",
                             "37-72 years" = "[37,73]"),
    
    mo_dipl = as.character(mo_dipl),                                            # niveau d'éducation de la mère
    mo_dipl_3cat = 
      fct_recode(mo_dipl,
                 "2 years or less after graduation" = "BEP/CAP/Highschool",
                 "2 years or less after graduation" = "1-2 years after graduation"), 
    mo_dipl_3cat = 
      fct_relevel(mo_dipl_3cat,
                  "≥5 years after graduation", 
                  "3-4 years after graduation", 
                  "2 years or less after graduation"), 
    mo_dipl_2cat = 
      fct_recode(mo_dipl_3cat,
                 "<5 years after graduation" = "3-4 years after graduation",
                 "<5 years after graduation" = "2 years or less after graduation"),
    fa_dipl_3cat =                                                              # niveau d'éducation du père
      fct_recode(fa_dipl,     
                 "2 years or less after graduation" = "1-2 years after graduation",
                 "2 years or less after graduation" = "BEP/CAP/Highschool"), 
    
    po_w_kg = po_w / 1000,                                                      # poids à la naissance
    po_w_kg_3cat = cut(po_w_kg,
                       include.lowest = TRUE,
                       right = FALSE,
                       dig.lab = 4,
                       breaks = c(0.9, 3, 3.5, 4.7)), 
    po_w_kg_3cat = fct_recode(po_w_kg_3cat,
                              "<3 Kg" = "[0.9,3)",
                              "3-3.4 Kg" = "[3,3.5)",
                              "≥3.5 Kg" = "[3.5,4.7]"),
    po_he_3cat = cut(po_he,                                                     # taille à la naissance
                     include.lowest = TRUE,
                     right = FALSE,
                     dig.lab = 4,
                     breaks = c(30, 50, 52, 60)),
    po_he_3cat = fct_recode(po_he_3cat,
                            "<50 cm" = "[30,50)",
                            "50-51 cm" = "[50,52)",
                            "≥52 cm" = "[52,60]"),
    
    ch_w_Y1_3cat = cut(ch_w_Y1,                                                 # poids à un an 
                       include.lowest = TRUE,
                       right = FALSE,
                       dig.lab = 4,
                       breaks = c(6.0, 8.5, 10, 14)),
    ch_w_Y1_3cat = fct_recode(ch_w_Y1_3cat,
                              "<8.5 Kg" = "[6,8.5)",
                              "8.5-9.9 Kg" = "[8.5,10)",
                              "≥10 Kg" = "[10,14]"),
    ch_he_Y1_3cat = cut(ch_he_Y1,                                               # taille à un an 
                        include.lowest = TRUE,
                        right = FALSE,
                        dig.lab = 4,
                        breaks = c(67.75, 75, 78, 83.75)),
    ch_he_Y1_3cat = fct_recode(ch_he_Y1_3cat,
                               "<75 cm" = "[67.75,75)",
                               "75-77.9 cm" = "[75,78)",
                               "≥78 cm" = "[78,83.75]"),
    
    ch_w_Y3_3cat = cut(ch_w_Y3,                                                 # poids à 3 ans
                       include.lowest = TRUE,
                       right = FALSE,
                       dig.lab = 2,
                       breaks = c(9, 13.85, 15.2333333333333, 19.6)), 
    ch_w_Y3_3cat = fct_recode(ch_w_Y3_3cat,
                              "9-13 Kg" = "[9,14)",
                              "14 Kg" = "[14,15)",
                              "15-20 Kg" = "[15,20]"),
    ch_he_Y3_3cat = cut(ch_he_Y3,                                               # taille à 3 ans
                        include.lowest = TRUE,
                        right = FALSE,
                        dig.lab = 2,
                        breaks = c(81, 93.75, 96.5, 106)),
    ch_he_Y3_3cat = fct_recode(ch_he_Y3_3cat, 
                               "81-93 cm" = "[81,94)",
                               "94-95 cm" = "[94,96)",
                               "96-106 cm" = "[96,1.1e+02]"), 
    
    mo_bmi_bepr_3cat = cut(mo_bmi_bepr,                                         # BMI pré grossesse de la mère
                           include.lowest = TRUE, 
                           right = FALSE,
                           dig.lab = 4,
                           breaks = c(16, 19, 24, 42)),
    mo_bmi_bepr_3cat = fct_recode(mo_bmi_bepr_3cat,
                                  "<19 Kg/m2" = "[16,19)",
                                  "19-23.9 Kg/m2" = "[19,24)",
                                  "≥24 Kg/m2" = "[24,42]"),
    
    po_gd_4cat = cut(po_gd,                                                     # durée de la gestation                          
                     include.lowest = TRUE,
                     right = FALSE,
                     dig.lab = 4,
                     breaks = c(28, 38, 40, 41, 42)), 
    po_gd_4cat = fct_recode(po_gd_4cat,                                     
                            "<38 weeks of amenorrhea" = "[28,38)",
                            "38-39 weeks of amenorrhea" = "[38,40)",
                            "40 weeks of amenorrhea" = "[40,41)",
                            ">40 weeks of amenorrhea" = "[41,42]"),
    
    ch_bf_duration_till48w_4cat = cut(ch_bf_duration_till48w,                   # durée de l'allaitement            
                                      include.lowest = TRUE,
                                      right = FALSE,
                                      dig.lab = 4,
                                      breaks = c(0, 1, 24, 47, 48)), 
    ch_bf_duration_till48w_4cat = fct_recode(ch_bf_duration_till48w_4cat,         
                                             "Not breastfed" = "[0,1)",
                                             "<24 weeks" = "[1,24)",
                                             "24-47 weeks" = "[24,47)",
                                             "Still breastfeed at 48 weeks" = "[47,48]"),   
  
    ch_antibio_Y1 = as.character(ch_antibio_Y1),                                # nombre de prescription d'antibiotiques
    ch_antibio_Y1_3cat = fct_recode(ch_antibio_Y1,
                                    "2 and more" = "2",
                                    "2 and more" = "3",
                                    "2 and more" = "4",
                                    "2 and more" = "5"),
    ch_antibio_Y1_2cat = fct_recode(ch_antibio_Y1,
                                    "No" = "0",
                                    "Yes" = "1",
                                    "Yes" = "2",
                                    "Yes" = "3",
                                    "Yes" = "4",
                                    "Yes" = "5"),
    
    ch_food_intro_Y1_3cat = fct_recode(ch_food_intro_Y1,                        # période de début de l'alimentation solide
                                       "Between 0 and 6 months old" = "Between 3 and 6 months old",
                                       "Between 0 and 6 months old" = "Between 0 and 3 months old"), 
   
    mo_alcool_preg_opt_3cat = fct_recode(mo_alcool_preg_opt,                    # consomation d'alcool pendant la grossesse
        "Did not drink" = "No drinking",
        "Drank only before finding out she was pregnant" = "<1drink/month when she didn’t know she was pregnant and did not drink when she knew she was pregnant (or did not answer)",
        "Drank only before finding out she was pregnant" = ">1drink/month when she didn’t know she was pregnant and did not drink when she knew she was pregnant",
        "Drank" = "<1drink/month when she knew she was pregnant",
        "Drank" = ">1drink/month when she knew she was pregnant"), 

    mo_alcool_preg_opt_4cat = fct_recode(
      mo_alcool_preg_opt,
      "Did not drink" = "No drinking",
      "Drank only before finding out she was pregnant" = "<1drink/month when she didn’t know she was pregnant and did not drink when she knew she was pregnant (or did not answer)",
      "Drank only before finding out she was pregnant" = ">1drink/month when she didn’t know she was pregnant and did not drink when she knew she was pregnant"), 
    
    ch_feces_age_w_Y1 = as.numeric(ch_feces_age_w_Y1),                          # age de l'enfant au prélévement de la selle                   
    ch_feces_age_w_Y1_4cat = cut(
      ch_feces_age_w_Y1,
      include.lowest = TRUE,
      right = FALSE,
      dig.lab = 4,
      breaks = c(42, 51, 53, 55, 71)),
    ch_feces_age_w_Y1_4cat = fct_recode(            
      ch_feces_age_w_Y1_4cat,
      "<51weeks" = "[42,51)",
      "51-52weeks" = "[51,53)",
      "53-54weeks" = "[53,55)",
      "≥55weeks" = "[55,71]"))


# Labellisation des variables ----
bdd = modify(bdd,{
  
  var_lab(ident) = "SEPAGES dentity"
  
  var_lab(mo_ethnicity) = "Maternal ethnicity"                                  # varibles parents
  var_lab(mo_ethnicity_2cat) = "Maternal ethnicity"
  var_lab(fa_ethnicity_2cat) = "Paternal ethinicity"
  
  var_lab(mo_dipl) = "Maternal education" 
  var_lab(mo_dipl_3cat) = "Maternal education"
  var_lab(mo_dipl_2cat) = "Maternal education"
  var_lab(fa_dipl) = "Paternal education"  
  var_lab(fa_dipl_3cat) = "Paternal education"  
  
  var_lab(mo_age) = "Maternal age at conception (years)"                               
  var_lab(mo_age_4cat) = "Maternal age at conception"
  var_lab(mo_age_3cat) = "Maternal age at conception"
  var_lab(fa_age) = "Paternal age at child birth (years)"
  var_lab(fa_age_3cat) = "Paternal age at child birth"
  
  var_lab(mo_coupl_pr) = "Woman living in couple during pregnancy"
  
  var_lab(mo_bmi_bepr) = "Maternal BMI before pregnancy (kg/m2)" 
  var_lab(mo_bmi_bepr_3cat) = "Maternal BMI before pregnancy" 
   
  var_lab(po_delmod) = "Delivery mode"                                          # variable accouchement          
  var_lab(mo_par) = "Maternal parity" 
  var_lab(mo_par_2cat) = "Maternal parity" 
  var_lab(po_gd) = "Gestational duration, completed weeks"    
  var_lab(po_gd_4cat) = "Gestational term, completed weeks" 
  
  var_lab(ch_sex) = "Child sex"                                                 # variable enfant
  
  var_lab(po_w) = "Birth weight (g)" 
  var_lab(po_w_kg) = "Birth weight (Kg)"
  var_lab(po_w_kg_3cat) = "Birth weight" 
  var_lab(po_he) = "Birth length (cm)" 
  var_lab(po_he_3cat) = "Birth length" 
  var_lab(ch_w_Y1) = "Weight at one year (Kg)" 
  var_lab(ch_he_Y1) = "Length at one year (cm)"
  var_lab(ch_w_Y1_3cat) = "Weight at one year"
  var_lab(ch_he_Y1_3cat) = "Length at one year"
  var_lab(ch_w_Y3) = "Weight at 3 years (Kg)"
  var_lab(ch_he_Y3) = "Length at 3 years (cm)"
  var_lab(ch_w_Y3_3cat) = "Weight at 3 years"
  var_lab(ch_he_Y3_3cat) = "Length at 3 years"
  
  var_lab(ch_bf_duration_till48w) = "Breastfeeding duration (weeks)"            ## nourriture et allaitement
  var_lab(ch_bf_duration_till48w_4cat) = "Breastfeeding duration"
  var_lab(ch_food_intro_Y1) = "Period of introduction of solid food 0-12 months old"
  var_lab(ch_food_intro_Y1_3cat) = "Period of introduction of solid food 0-12 months old"
  
  var_lab(mo_pets) = "Presence of pets"                                         ## pets
  
  var_lab(ch_care_main_6m_12m_opt2_2c) = "Main mode of child care up to 1 year" ## mode of care
  var_lab(ch_care_main_6m_12m_opt2_3c) = "Main mode of child care up to 1 year"
  
  var_lab(ch_antibio_Y1) = "Number of antibiotics use between 0-12 months old"  ## expo
  var_lab(ch_antibio_Y1_3cat) = "Number of antibiotics use between 0-12 months old"
  var_lab(ch_antibio_Y1_2cat) = "Antibiotics use between 0-12 months old"
  
  var_lab(mo_tob_gr_anyt_yn_n2) = "Maternal active smoking during pregnancy"    
  var_lab(Mo_ETS_anyT_yn1_opt) = "Maternal passive smoking during pregnancy"
  var_lab(ch_ETS_12m_opt36m) = "Child passive smoking during the first year of life"
  var_lab(ch_tabacco_passive_up_to_Y1) = "Child perinatal passive smoking up to one year"
  var_lab(ch_tabacco_total_Y1) = "Child perinatal exposure to tabacco up to one year"
  
  var_lab(mo_alcool_preg) = "Maternal alcohol consumption during pregnancy"     
  var_lab(mo_alcool_preg_opt) = "Maternal alcohol consumption during pregnancy"
  var_lab(mo_alcool_preg_opt_3cat) = "Maternal alcohol consumption during pregnancy"
  var_lab(mo_alcool_preg_opt_4cat) = "Maternal alcohol consumption during pregnancy"
  
  var_lab(mo_hadanxscore_grt3_imp) = "Maternal anxiety score during the third trimester of pregnancy"
  var_lab(mo_haddepscore_grt3_imp) = "Maternal depression score during the third trimester of pregnancy"
  var_lab(mo_hadtotscore_grt3_imp) = "Maternal anxiety and depression score during the third trimester of pregnancy"
  var_lab(mo_hadanxscore_y1_imp) =  "Maternal anxiety score when child is one year old" 
  var_lab(mo_haddepscore_y1_imp) = "Maternal depression score when child is one year old" 
  var_lab(mo_hadtotscore_y1_imp) = "Maternal anxiety and depression score when child is one year old"
  
  var_lab(home_learning_y3) = "Learning materials score at 3 years of age (HOME questionnaire)"
  var_lab(home_language_y3) = "Language stimulation score at 3 years of age (HOME questionnaire)"
  var_lab(home_academic_y3) = "Academic stimulation score at 3 years of age (HOME questionnaire)"
  var_lab(home_variety_y3) = "Variety score at 3 years of age (HOME questionnaire)"
  var_lab(home_total_y3) = "HOME questionnaire total score"

  var_lab(ch_feces_age_w_Y1) = "Child age at stool collection (weeks)"          # variables microbiote    
  var_lab(ch_feces_age_w_Y1_4cat) = "Child age at stool collection"
  var_lab(ch_feces_RUN_Y1) = "MiSeq analytic batch effect"
  
  var_lab(ch_cbclintscore_y2) = "Child interne score at 2 years (CBCL)"         # variabls neuro
  var_lab(ch_cbclextscore_y2) = "Child externe score at 2 years (CBCL)"
  
  var_lab(ch_socawar_y3) = "Child social awareness score at 3 years (SRS)"    
  var_lab(ch_soccog_y3) = "Child social cognition score at 3 years (SRS)" 
  var_lab(ch_soccom_y3) = "Child social communication score at 3 years (SRS)"   
  var_lab(ch_socmot_y3) = "Child social motivation score at 3 years (SRS)"  
  var_lab(ch_RRB_y3) = "Child restricted interests and repetitive behavior score at 3 years (SRS)"  
  var_lab(ch_SRStotal_y3) = "Child SRS total score at 3 years (SRS)"    
  
  var_lab(ch_briefpinhibit_y3) = "Child inhibition score at 3 years (BRIEF-P)"  
  var_lab(ch_briefpshift_y3) = "Child shift score at 3 years (BRIEF-P)"                               
  var_lab(ch_briefpemocontrol_y3) = "Child emotional control score at 3 years (BRIEF-P)"  
  var_lab(ch_briefpworkmemo_y3) = "Child working memory score at 3 years (BRIEF-P)"  
  var_lab(ch_briefpplan_y3) = "Child plan/organization score at 3 years (BRIEF-P)"  
  
  var_lab(ch_verbal_comprehension_IQ_Y3) = "Verbal comprehension score at 3 years (IQ test)"
  var_lab(ch_visuospatiale_IQ_Y3) = "Visuospatiale score at 3 years (IQ test)"
  var_lab(ch_work_memory_IQ_Y3) = "Work memory score at 3 years (IQ test)"
  var_lab(ch_total_IQ_Y3) = "Total IQ score at 3 years"
  
  var_lab(ch_age_CBCL_Y2) = "Child age at the CBCL 2 years evaluation"
  var_lab(ch_age_Vineland_Y2) = "Child age at the Vineland 2 years evaluation (cy2hap2)"
  var_lab(ch_age_Vineland_Y2_bis) = "Child age at the Vineland 2 years evaluation (cy2hap3)"
  var_lab(ch_age_IQ_Y3) = "Child age at the IQ 3 years evaluation"
  var_lab(ch_age_SRS_BRIEFP_Y3) = "Child age at the SRS 3 years evaluation"
  
    })

# Export ----
bdd$ident <- as.character(bdd$ident)
bdd <- left_join(bdd, bdd_genera_imp, by ="ident")

save.image(
  file = "1_intermediate_data/1_data_cleaning_AD_gumme.RData")
