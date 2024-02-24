# Bazar
# Aline 
# 18/12/2023


# Vérification que les fonctions ont fonctionné correctement ----
test <- with(data = bdd_final_imp, 
             exp = lm(ch_cbclintscore_y2 ~ 
                        ch_feces_rel_p4_Y1 +
                        
                        ch_age_CBCL_Y2 +
                        
                        po_w_kg_3cat + 
                        po_he_3cat + 
                        mo_dipl_2cat +
                        mo_age +
                        mo_bmi_bepr_3cat +    
                        ch_sex +
                        mo_par_2cat + 
                        ch_bf_duration_till48w_4cat +
                        po_gd +
                        po_delmod + 
                        ch_food_intro_Y1_3cat +
                        mo_pets +
                        ch_antibio_Y1_2cat +
                        #home_total_y3 +                                             # on ne met pas home Y3 dans le modèle Y2
                        mo_hadtotscore_grt3_imp +   
                        mo_tob_gr_anyt_yn_n2+
                        ch_tabacco_passive_up_to_Y1 +
                        ch_care_main_12m_opt2_2c))
test_pooled <- pool(test)
#summary(test)
summary(test_pooled)
test %>% tbl_regression() %>% plot()


phylum_label <- bdd %>% select(contains("ch_feces_rel_p")) %>% var_label() %>% unlist()
tbl_merge(
  tbls = 
    list(
      tbl_1 = 
        bdd %>% 
        filter(!is.na(ch_feces_rel_g1_Y1)) %>%
        select(contains("ch_feces_rel_p")) %>% 
        tbl_summary(
          type = list(everything() ~ "continuous"), 
          statistic = list(everything() ~ "{median} ({p25}, {p75})"), 
          digits = list(all_continuous() ~ c(2, 1, 1))), 
      tbl_2 = 
        bdd %>%
        filter(!is.na(ch_feces_rel_g1_Y1)) %>%
        select(contains("ch_feces_rel_p")) %>% 
        mutate_all(~ ifelse(.>0, "Yes", "No")) %>%
        set_label(phylum_label) %>%
        tbl_summary(
          type = list(everything() ~ "categorical"))), 
  tab_spanner = c("**Continuous**", "**Categorical (Y/N)**"))
rm(phylum_label)

tbl_merge(
  tbls = 
    list(
      tbl_1 = 
        bdd %>% 
        filter(!is.na(ch_feces_rel_g1_Y1)) %>%
        select(contains("ch_feces_rel_c")) %>% 
        tbl_summary(
          type = list(everything() ~ "continuous"), 
          statistic = list(everything() ~ "{median} ({p25}, {p75})"), 
          digits = list(all_continuous() ~ c(2, 1, 1))), 
      tbl_2 = 
        bdd %>%
        filter(!is.na(ch_feces_rel_g1_Y1)) %>%
        select(contains("ch_feces_rel_c")) %>% 
        mutate_all(~ ifelse(.>0, "Yes", "No")) %>%
        set_label(class_label) %>%
        tbl_summary(
          type = list(everything() ~ "categorical"))), 
  tab_spanner = c("**Continuous**", "**Categorical (Y/N)**"))


# table taxonomique ----
test <- taxa_table_Y1 %>% 
  select(-ch_feces_ASV_ID_Y1) %>% 
  unique() %>%
  mutate(
    ch_feces_phylum_ASVbased_Y1 = fct_relevel(
      ch_feces_phylum_ASVbased_Y1,
      "Actinobacteria", "Firmicutes", "Bacteroidetes", "Proteobacteria",
      "Verrucomicrobia", "Fusobacteria", "Tenericutes", "Candidatus_Saccharibacteria",
      "Cyanobacteria_Chloroplast", "Euryarchaeota")) %>%
  arrange(ch_feces_phylum_ASVbased_Y1)


# genres d'interet ----
bdd %>%
  select(contains("ch_feces_rel_g")) %>%
  select(!contains("_ln")) %>%
  var_label()

ch_feces_rel_g5_Y1   # Akkermansia
ch_feces_rel_g15_Y1  # Enterococus
ch_feces_rel_g33_Y1  # Lactobacillus
ch_feces_rel_g34_Y1  # Klebsiella
ch_feces_rel_g138_Y1 # Campylobacter
# Faacalibacterium

bdd %>%
  filter(!is.na(ch_feces_rel_g1_Y1)) %>%
  select(ch_feces_rel_g5_Y1,   # Akkermansia
         ch_feces_rel_g15_Y1,  # enterococus
         ch_feces_rel_g33_Y1,  # Lactobacillus
         ch_feces_rel_g34_Y1,  # Klebsiella
         ch_feces_rel_g138_Y1) %>%  # Campylobacter
  tbl_summary()


verif_distrib(data = bdd, var = ch_feces_rel_g5_Y1)
verif_distrib(data = bdd, var = ch_feces_rel_g15_Y1)
verif_distrib(data = bdd, var = ch_feces_rel_g33_Y1)
verif_distrib(data = bdd, var = ch_feces_rel_g34_Y1)
verif_distrib(data = bdd, var = ch_feces_rel_g138_Y1)

bdd %>% 
  mutate(
    ch_feces_rel_g5_Y1_cat = ifelse(ch_feces_rel_g5_Y1 > 0, "Yes", "No"), 
    ch_feces_rel_g15_Y1_cat = ifelse(ch_feces_rel_g15_Y1 > 0, "Yes", "No"), 
    ch_feces_rel_g33_Y1_cat = ifelse(ch_feces_rel_g33_Y1 > 0, "Yes", "No"), 
    ch_feces_rel_g34_Y1_cat = ifelse(ch_feces_rel_g34_Y1 > 0, "Yes", "No"), 
    ch_feces_rel_g138_Y1_cat = ifelse(ch_feces_rel_g138_Y1 > 0, "Yes", "No")) %>%
  filter(!is.na(ch_feces_rel_g1_Y1)) %>%
  select(ch_feces_rel_g5_Y1, ch_feces_rel_g5_Y1_cat, 
         ch_feces_rel_g15_Y1, ch_feces_rel_g15_Y1_cat, 
         ch_feces_rel_g33_Y1, ch_feces_rel_g33_Y1_cat, 
         ch_feces_rel_g34_Y1, ch_feces_rel_g34_Y1_cat, 
         ch_feces_rel_g138_Y1, ch_feces_rel_g138_Y1_cat) %>%
  tbl_summary(
    type = list(all_categorical() ~ "categorical"))



# Création de variables genres catégorielles ----
tbl_1 <- bdd %>% 
  filter(!is.na(ch_feces_rel_g1_Y1)) %>%
  select(ident, contains("ch_feces_rel_g")) %>%
  select(!contains("_ln")) %>%
  select(ident, order(-sapply(., median))) %>%
  tbl_summary(
    include = -ident,
    type = list(everything() ~ "continuous"), 
    statistic = list(everything() ~ "{median} ({p25}, {p75})"), 
    digits = list(all_continuous() ~ c(2, 1, 1)))

genera_labels <- 
  bdd %>%
  select(ident, 
         contains("ch_feces_rel_g"))%>% 
  select(!contains("_ln")) %>% 
  var_label() %>% 
  unlist()

bdd_genera <- bdd %>%
  select(ident, 
         contains("ch_feces_rel_g"))%>% 
  select(!contains("_ln")) %>%
  set_label(genera_labels)

bdd_genera_2 <- bdd %>%
  select(ident, 
         contains("ch_feces_rel_g"))%>% 
  select(!contains("_ln")) %>%
  mutate_at(vars(-ident), ~ ifelse(.>0, "Yes", "No")) %>%
  #rename_at(-1, ~ paste0(., "_cat")) %>%
  set_label(genera_labels)

tbl_2 <- bdd_genera_2 %>%
  select(-ident) %>%
  filter(!is.na(ch_feces_rel_g1_Y1)) %>%
  #select(contains("_cat")) %>%
  tbl_summary(
    type = list(everything() ~ "categorical"))

bdd_genera_2 <- bdd_genera_2 %>% rename_at(-1, ~ paste0(., "_cat")) 


bdd_genera <- left_join(bdd_genera, bdd_genera_2, by ="ident")

tbl_descrip_genera <- tbl_merge(
  tbls = list(tbl_1, tbl_2), 
  tab_spanner = c("**Continuous**", "**Categorical (Y/N)**"))
tbl_descrip_genera

rm(bdd_genera_2, genera_labels, tbl_1, tbl_2)



# Test taxonomic packages
## LEFSE ----
BiocManager::install("lefser")
library(lefser)


data(zeller14)
zeller14 <- zeller14[, zeller14$study_condition != "adenoma"]
table(zeller14$study_condition)
table(zeller14$age_category)
table(zeller14$age_category, zeller14$study_condition)
res <- lefser(zeller14, groupCol = "study_condition", blockCol = "age_category")
lefserPlot(res)


## ANCOM ----
BiocManager::install("ANCOMBC")
library(ANCOMBC)


data(atlas1006, package = "microbiome")
tse = mia::makeTreeSummarizedExperimentFromPhyloseq(atlas1006)

# subset to baseline
tse = tse[, tse$time == 0]

# Re-code the bmi group
tse$bmi = recode(tse$bmi_group,
                 obese = "obese",
                 severeobese = "obese",
                 morbidobese = "obese")
# Subset to lean, overweight, and obese subjects
tse = tse[, tse$bmi %in% c("lean", "overweight", "obese")]

# Note that by default, levels of a categorical variable in R are sorted 
# alphabetically. In this case, the reference level for `bmi` will be 
# `lean`. To manually change the reference level, for instance, setting `obese`
# as the reference level, use:
tse$bmi = factor(tse$bmi, levels = c("obese", "overweight", "lean"))
# You can verify the change by checking:
# levels(sample_data(tse)$bmi)

# Create the region variable
tse$region = recode(as.character(tse$nationality),
                    Scandinavia = "NE", UKIE = "NE", SouthEurope = "SE", 
                    CentralEurope = "CE", EasternEurope = "EE",
                    .missing = "unknown")

# Discard "EE" as it contains only 1 subject
# Discard subjects with missing values of region
tse = tse[, ! tse$region %in% c("EE", "unknown")]

print(tse)