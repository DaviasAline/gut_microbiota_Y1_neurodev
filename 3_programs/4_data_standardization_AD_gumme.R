## Aline Davias
## Standardisation des variables microbiote d'intéret
## 27/03/2024


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
source("3_programs/5_functions_AD_gumme.R")

# Création des vecteurs de variables ----
protocol_vars <- c("ch_feces_age_w_Y1", "ch_feces_RUN_Y1")
# covariates_gm <- c("po_w_kg_3cat", 
#                    "po_he_3cat", 
#                    "mo_dipl_2cat", 
#                    "mo_age", 
#                    "mo_bmi_bepr_3cat", 
#                    "ch_sex", 
#                    "mo_par_2cat", 
#                    "ch_bf_duration_till48w_4cat", 
#                    "po_gd", 
#                    "po_delmod", 
#                    "ch_food_intro_Y1_3cat", 
#                    "mo_pets", 
#                    "ch_antibio_Y1_2cat", 
#                    "mo_hadtotscore_grt3_imp",     
#                    "mo_tob_gr_anyt_yn_n2", 
#                    "ch_tabacco_passive_up_to_Y1",
#                    "ch_care_main_12m_opt2_2c")

# Standardisation ----
bdd_a_std <- bdd_final_imp_1 %>%
  select(ident, 
         all_of(microbiote_vec), 
         all_of(protocol_vars))

bdd_a_std <- bdd_a_std %>% mutate(across(c(microbiote_vec), ~as.numeric(unclass(.))))
bdd_a_std_long <- bdd_a_std %>% 
  pivot_longer(cols = microbiote_vec) %>%
  select(ident, name, value, everything()) %>% 
  drop_all_labels()


bdd_std_long <- bdd_a_std_long %>%                    #' The function is designed to be called within a `mutate` call on grouped data where each group represents a different exposure.
  mutate(
    name = as.character(name)) %>%
  mutate(
    var_stad = standardise(                   
      var_to_std = "value",                #' @param var_to_std A character string with the name of the variable to be standardised.
      protocol_vars = protocol_vars,     #' @param protocol_vars A character vector of names of potential protocol variables.
      folder = "4_output/std/",
    ),
    .by = name
  )


bdd_std_long <- bdd_std_long %>%
  rename(microbiote = name) %>%
  pivot_longer(cols = -c("ident", "microbiote", "ch_feces_age_w_Y1", "ch_feces_RUN_Y1")) %>%
  rename(standardization = name) %>%
  mutate(standardization = fct_recode(standardization, 
                                      "Raw" = "value",
                                      "Standardized" = "var_stad"))
rm(bdd_a_std, bdd_a_std_long)
  

# Visu ----
bdd_std_long %>%
  filter(microbiote %in% c("ch_feces_SpecRich_5000_ASV_Y1_10", 
                           "ch_feces_Shannon_5000_ASV_Y1", 
                           "ch_feces_SpecRich_10000_ASV_Y1_10", 
                           "ch_feces_Shannon_10000_ASV_Y1", 
                           "ch_feces_rel_p1_Y1_10",
                           "ch_feces_rel_p2_Y1_10",
                           "ch_feces_rel_p3_Y1_10",
                           "ch_feces_rel_p4_Y1_10" )) %>%
  ggplot() +
  aes(x = value, y = standardization, fill = ch_feces_RUN_Y1) +
  geom_boxplot() +
  scale_fill_hue(direction = 1) +
  coord_flip() +
  theme_minimal() +
  facet_wrap(~microbiote, ncol = 4, scales = "free_y") +
  theme_lucid()


# pivot_wider ----
## Récup variables brutes ----
bdd_std_large_raw <- bdd_std_long %>%
  filter(standardization == "Raw") %>%
  pivot_wider(names_from = standardization, 
              values_from = value) 

bdd_std_large_raw <- bdd_std_large_raw %>%
  pivot_wider(names_from = microbiote, values_from = Raw)

## Récup variables std ----
bdd_std_large_std <- bdd_std_long %>%
  filter(standardization == "Standardized") %>%
  pivot_wider(names_from = standardization, 
              values_from = value) 

bdd_std_large_std <- bdd_std_large_std %>%
  pivot_wider(names_from = microbiote, values_from = Standardized) 

colnames(bdd_std_large_std) <- colnames(bdd_std_large_std) %>%
  str_replace_all("_imp_log_Y1", "_imp_log_std_Y1") %>%
  str_replace_all("_Y1_10", "_std_Y1_10") %>%
  str_replace_all("ASV_Y1", "ASV_std_Y1")

microbiote_std_vec <- colnames(select(bdd_std_large_std, -"ident", -"ch_feces_age_w_Y1", -"ch_feces_RUN_Y1"))

bdd_std_large <- left_join(bdd_std_large_raw, 
                           bdd_std_large_std, 
                           by = c("ident", "ch_feces_age_w_Y1", "ch_feces_RUN_Y1"))
rm(bdd_std_large_raw, bdd_std_large_std)



## Visu corrélation 
bdd_std_long <- bdd_std_long %>%
  pivot_wider(names_from = standardization, values_from = value)

resultats_correlation <- bdd_std_long %>% 
  group_by(microbiote) %>% 
  summarise(correlation = cor(Raw, Standardized, use = "complete.obs")) %>%
  ungroup() 

resultats_correlation %>% View()

# Merge avec la base de données principales ----
bdd_final_imp_1$ident <- as.character(bdd_final_imp_1$ident)
bdd_final_imp_1 <- bdd_final_imp_1 %>% drop_all_labels()
bdd_final_imp_1 <- left_join(bdd_std_large, bdd_final_imp_1, 
                             by = c("ident", protocol_vars, microbiote_vec))

bdd_final_imp_1_sensi_seuil$ident <- as.character(bdd_final_imp_1_sensi_seuil$ident)
bdd_final_imp_1_sensi_seuil <- bdd_final_imp_1_sensi_seuil %>% drop_all_labels()
bdd_final_imp_1_sensi_seuil <- left_join(bdd_std_large, bdd_final_imp_1_sensi_seuil, 
                                         by = c("ident", protocol_vars, microbiote_vec))



# Export ----
save(list = ls(),
     file = "1_intermediate_data/4_data_standardization_AD_gumme.RData")
