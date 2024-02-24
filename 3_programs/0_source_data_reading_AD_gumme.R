# 0_source_data_reading_AD_gumme
# Aline Davias
# 25/08/2023

# Chargement des packages ----
library(readxl)
library(haven)
library(tidyverse)
library(expss)

# Chargement des données ----
metadata <- read_sas("C:/Users/Aline/OneDrive - etu.univ-grenoble-alpes.fr/Documents/5. R projects/gut_microbiota_Y1_neurodev/0_source_data/base_aline_230925.sas7bdat", NULL)                                   
alphadiv_Y1 <-  read_labelled_csv("C:/Users/Aline/OneDrive - etu.univ-grenoble-alpes.fr/Documents/5. R projects/gut_microbiota_Y1_neurodev/0_source_data/alpha_diversity_ASVbased_Y1_labelled_AD_20220504_26.csv")
microbiote_Y1 <- read_labelled_csv("C:/Users/Aline/OneDrive - etu.univ-grenoble-alpes.fr/Documents/5. R projects/gut_microbiota_Y1_neurodev/0_source_data/gut_microbiota_ASVbased_Y1_labelled_AD_20220504_7239.csv")
age_feces_Y1 <- read_labelled_csv("C:/Users/Aline/OneDrive - etu.univ-grenoble-alpes.fr/Documents/5. R projects/gut_microbiota_Y1_neurodev/0_source_data/age_feces_collection_Y1_labelled_AD_20220504_4.csv")
atb_Y1 <- read_labelled_csv("C:/Users/Aline/OneDrive - etu.univ-grenoble-alpes.fr/Documents/5. R projects/gut_microbiota_Y1_neurodev/0_source_data/antibiotics_Y1_labelled_AD_20220314_9.csv")
solidfood_Y1 <- read_labelled_csv("C:/Users/Aline/OneDrive - etu.univ-grenoble-alpes.fr/Documents/5. R projects/gut_microbiota_Y1_neurodev/0_source_data/solidfood_Y1_labelled_AD_20220302_7.csv")
weight_length_Y1 <- read_labelled_csv("C:/Users/Aline/OneDrive - etu.univ-grenoble-alpes.fr/Documents/5. R projects/gut_microbiota_Y1_neurodev/0_source_data/weight_length_Y1_labelled_AD_20220301_2.csv")

taxa_table_Y1 <- read_csv("C:/Users/Aline/OneDrive - etu.univ-grenoble-alpes.fr/Documents/5. R projects/gut_microbiota_Y1_neurodev/0_source_data/taxa_table_ASVbased_Y1_AD_20220504_8.csv")  # ne pas inclure, table de correspondance taxonomique
load("C:/Users/Aline/OneDrive - etu.univ-grenoble-alpes.fr/Documents/5. R projects/gut_microbiota_Y1_neurodev/0_source_data/ch_mode_care_AD_20231005_49.RData") 
bdd_care <- bdd_care %>% select(ident, contains("main"))

var_lab(metadata$ident) <- NULL
var_lab(alphadiv_Y1$ident) <- NULL
var_lab(microbiote_Y1$ident) <- NULL
var_lab(age_feces_Y1$ident) <- NULL
var_lab(atb_Y1$ident) <- NULL
var_lab(solidfood_Y1$ident) <- NULL
var_lab(weight_length_Y1$ident) <- NULL
var_lab(bdd_care$ident) <- NULL

# Regroupement des données ----
bdd <- 
  reduce(
    list(metadata, alphadiv_Y1, microbiote_Y1, age_feces_Y1, atb_Y1, solidfood_Y1, weight_length_Y1, bdd_care), 
    full_join, 
    by = "ident")

save(bdd, 
     file = "C:/Users/Aline/OneDrive - etu.univ-grenoble-alpes.fr/Documents/5. R projects/gut_microbiota_Y1_neurodev/1_intermediate_data/0_source_data_reading_AD_gumme.RData")
