
library(broom)         
library(effectsize)   
library(dplyr)          

explanatory_names <- c("Specific richness", "Shannon diversity", 
                       "Firmicutes", "Actinobacteria", "Bacteroidetes", "Proteobacteria", 
                       genera_var_labels)


results <- data.frame()

for (outcome in outcomes) {
  
  if (outcome %in% c("ch_cbclintscore_y2", "ch_cbclextscore_y2")) {
    covariate_specific <- covariates_CBCL
  } else if (outcome %in% c("ch_SRStotal_y3", "ch_briefpinhibit_y3", "ch_briefpshift_y3", "ch_briefpemocontrol_y3", "ch_briefpworkmemo_y3", "ch_briefpplan_y3")) {
    covariate_specific <- covariates_SRS_BRIEF
  } else {
    covariate_specific <- covariates_IQ
  }
  
  for (i in seq_along(explanatory)) {
    
    explanatory_var <- explanatory[i]
    explanatory_name <- explanatory_names[i]
    
    formula_complet <- as.formula(paste(outcome, "~", explanatory_var, "+", paste(covariate_specific, collapse = " + ")))
    formula_reduit <- as.formula(paste(outcome, "~", paste(covariate_specific, collapse = " + ")))
    
    model_complet <- lm(formula_complet, data = bdd_final_imp_1)
    model_reduit <- lm(formula_reduit, data = bdd_final_imp_1)
    
    model_summary <- tidy(model_complet, conf.int = TRUE)
    
    model_summary_exp <- model_summary %>%
      filter(term == explanatory_var) %>%
      select(term, estimate, conf.low, conf.high, p.value) %>%
      mutate(outcome = outcome)
    
    R2_complet <- summary(model_complet)$r.squared
    R2_reduit <- summary(model_reduit)$r.squared
    
    Cohens_f2 <-  (R2_complet - R2_reduit) / (1 - R2_complet)

    model_summary_exp <- model_summary_exp %>%
      mutate(Cohens_f2, term_2 = explanatory_name)
    
    results <- rbind(results, model_summary_exp)
  }
}
rm(outcome, covariate_specific, explanatory_var, explanatory_name, model_complet, model_reduit, model_summary, f2_value, model_summary_exp, i)

results <- results %>%
  mutate(
    Beta = format(round(estimate, 2), digits = 2),
    conf.low = format(round(conf.low, 2), digits = 2),
    conf.high = format(round(conf.high, 2), digits = 2),
    "95% CI" = paste(conf.low, conf.high, sep = ","), 
    p.value = case_when(p.value < 0.001 ~ format(round(p.value, 4), nsmall = 4),
                        p.value < 0.01 ~ format(round(p.value, 3), nsmall = 3),
                        p.value > 0.01 ~ format(round(p.value, 2), nsmall = 2)), 
    Cohens_f2 = case_when(Cohens_f2 <0 ~ "<0",
                          Cohens_f2> 0 & Cohens_f2 < 0.001 ~ "<0.001",
                          Cohens_f2> 0 & Cohens_f2 < 0.01 ~ format(round(Cohens_f2, 3), nsmall = 3),
                          Cohens_f2> 0 & Cohens_f2 > 0.01 ~ format(round(Cohens_f2, 2), nsmall = 2), 
                          .default = as.character(Cohens_f2))) %>%
  select(outcome, term, term_2, Beta, `95% CI`, `p.value`, Cohens_f2)

table_2_bis <- results %>%
  filter(term_2 %in% c("Specific richness", "Shannon diversity")) %>%
  rename(Explanatory = term_2) %>%
  select(-term)

table_4_bis <- results %>%
  filter(term_2 %in% c("Firmicutes", "Actinobacteria", "Bacteroidetes", "Proteobacteria")) %>%
  rename(Explanatory = term_2) %>%
  select(-term)

table_S4_bis <- results %>%
  filter(!term_2 %in% c("Specific richness", "Shannon diversity", 
                        "Firmicutes", "Actinobacteria", 
                        "Bacteroidetes", "Proteobacteria")) %>%
  select(-term) %>%
  pivot_wider(names_from = outcome, 
              values_from = c("Beta", "95% CI", "p.value", "Cohens_f2")) %>%
  mutate(
    term_2 = fct_recode(term_2, 
                             "Clostridium IV" = "Clostridium_IV",
                             "Clostridium sensu stricto" = "Clostridium_sensu_stricto",
                             "Clostridium XlVa" = "Clostridium_XlVa",
                             "Clostridium XVIII" = "Clostridium_XVIII",
                             "Erysipelotrichaceae incertae sedis" = "Erysipelotrichaceae_incertae_sedis",
                             "Escherichia and Shigella" = "Escherichia_Shigella",
                             "Lachnospiracea incertae sedis" = "Lachnospiracea_incertae_sedis",
                             "Ruminococcus 2" = "Ruminococcus2",
                             "Saccharibacteria genera incertae sedis" = "Saccharibacteria_genera_incertae_sedis"))%>%
  select("Gut microbiota parameters" = term_2, 
         contains("ch_cbclintscore_y2"), 
         contains("ch_cbclextscore_y2"),
         contains("ch_SRStotal_y3"), 
         contains("ch_briefpinhibit_y3"),
         contains("ch_briefpshift_y3"), 
         contains("ch_briefpemocontrol_y3"), 
         contains("ch_briefpworkmemo_y3"), 
         contains("ch_briefpplan_y3"), 
         contains("ch_WPPSI_verbal_comprehension_cor_Y3"), 
         contains("ch_WPPSI_visuospatiale_cor_Y3"), 
         contains("ch_WPPSI_work_memory_cor_Y3"), 
         contains("ch_WPPSI_total_cor_Y3"))

write.xlsx(table_2_bis, file = "4_output/table_2_bis.xlsx")
write.xlsx(table_4_bis, file = "4_output/table_4_bis.xlsx")
write.xlsx(table_S4_bis, file = "4_output/table_S4_bis.xlsx")


table_2_bis %>% filter(p.value<0.05 & Cohens_f2>0.02)
table_4_bis %>% filter(p.value<0.05 & Cohens_f2>0.02)
table_S4_bis %>% filter(p.value<0.05 & Cohens_f2>0.02)
results %>% filter(p.value<0.05 & Cohens_f2>0.019) %>% View()
results %>% filter(p.value<0.05 & Cohens_f2>0.12) %>% View()



## ch_cbclintscore_y2 ----
results %>% filter(p.value<0.05) %>% filter(outcome == "ch_cbclintscore_y2")%>% select(term)
model_reduit_cbcl_int <- lm(
  ch_cbclintscore_y2 ~
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
    home_total_y3 + 
    mo_hadtotscore_grt3_imp + 
    mo_tob_gr_anyt_yn_n2 +
    ch_tabacco_passive_up_to_Y1 + 
    ch_care_main_12m_opt2_2c,
  data = bdd_final_imp_1)

model_complet_cbcl_int <- lm(
  ch_cbclintscore_y2 ~
    ch_feces_rel_g19_imp_log_std_Y1 +
    ch_feces_rel_g21_imp_log_std_Y1 +
    ch_feces_rel_g36_imp_log_std_Y1 +
    ch_feces_rel_g45_imp_log_std_Y1 +
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
    home_total_y3 + 
    mo_hadtotscore_grt3_imp + 
    mo_tob_gr_anyt_yn_n2 +
    ch_tabacco_passive_up_to_Y1 + 
    ch_care_main_12m_opt2_2c,
  data = bdd_final_imp_1)

## ch_cbclextscore_y2 ----
results %>% filter(p.value<0.05) %>% filter(outcome == "ch_cbclextscore_y2") %>% select(term)

model_reduit_cbcl_ext <- lm(
  ch_cbclextscore_y2 ~
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
    home_total_y3 + 
    mo_hadtotscore_grt3_imp + 
    mo_tob_gr_anyt_yn_n2 +
    ch_tabacco_passive_up_to_Y1 + 
    ch_care_main_12m_opt2_2c,
  data = bdd_final_imp_1)

model_complet_cbcl_ext <- lm(
  ch_cbclextscore_y2 ~
    ch_feces_rel_p4_std_Y1_10 +
    ch_feces_rel_g1_imp_log_std_Y1 +
    ch_feces_rel_g13_imp_log_std_Y1 +
    ch_feces_rel_g15_imp_log_std_Y1 +
    ch_feces_rel_g19_imp_log_std_Y1 +
    ch_feces_rel_g34_imp_log_std_Y1 +
    ch_feces_rel_g39_imp_log_std_Y1 +
    ch_feces_rel_g43_imp_log_std_Y1 +
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
    home_total_y3 + 
    mo_hadtotscore_grt3_imp + 
    mo_tob_gr_anyt_yn_n2 +
    ch_tabacco_passive_up_to_Y1 + 
    ch_care_main_12m_opt2_2c,
  data = bdd_final_imp_1)



## ch_SRStotal_y3 ----
results %>% filter(p.value<0.05) %>% filter(outcome == "ch_SRStotal_y3") %>% select(term)

model_reduit_srs <- lm(
  ch_SRStotal_y3 ~
    ch_age_SRS_BRIEFP_Y3 + 
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
    home_total_y3 + 
    mo_hadtotscore_grt3_imp + 
    mo_tob_gr_anyt_yn_n2 +
    ch_tabacco_passive_up_to_Y1 + 
    ch_care_main_12m_opt2_2c,
  data = bdd_final_imp_1)

model_complet_srs <- lm(
  ch_SRStotal_y3 ~
    ch_feces_rel_p4_std_Y1_10 +
    ch_age_SRS_BRIEFP_Y3 + 
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
    home_total_y3 + 
    mo_hadtotscore_grt3_imp + 
    mo_tob_gr_anyt_yn_n2 +
    ch_tabacco_passive_up_to_Y1 + 
    ch_care_main_12m_opt2_2c,
  data = bdd_final_imp_1)

## ch_briefpinhibit_y3 ----
results %>% filter(p.value<0.05) %>% filter(outcome == "ch_briefpinhibit_y3") %>% select(term)

model_reduit_briefpinhibit <- lm(
  ch_briefpinhibit_y3 ~
    ch_age_SRS_BRIEFP_Y3 + 
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
    home_total_y3 + 
    mo_hadtotscore_grt3_imp + 
    mo_tob_gr_anyt_yn_n2 +
    ch_tabacco_passive_up_to_Y1 + 
    ch_care_main_12m_opt2_2c,
  data = bdd_final_imp_1)

model_complet_briefpinhibit <- lm(
  ch_briefpinhibit_y3 ~
    ch_feces_rel_p4_std_Y1_10 +
    ch_age_SRS_BRIEFP_Y3 + 
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
    home_total_y3 + 
    mo_hadtotscore_grt3_imp + 
    mo_tob_gr_anyt_yn_n2 +
    ch_tabacco_passive_up_to_Y1 + 
    ch_care_main_12m_opt2_2c,
  data = bdd_final_imp_1)

## ch_briefpshift_y3 ----
results %>% filter(p.value<0.05) %>% filter(outcome == "ch_briefpshift_y3") %>% select(term)

model_reduit_briefpshift <- lm(
  ch_briefpshift_y3 ~
    ch_age_SRS_BRIEFP_Y3 + 
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
    home_total_y3 + 
    mo_hadtotscore_grt3_imp + 
    mo_tob_gr_anyt_yn_n2 +
    ch_tabacco_passive_up_to_Y1 + 
    ch_care_main_12m_opt2_2c,
  data = bdd_final_imp_1)

model_complet_briefpshift <- lm(
  ch_briefpshift_y3 ~
    ch_feces_rel_g47_imp_log_std_Y1 +
    ch_age_SRS_BRIEFP_Y3 + 
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
    home_total_y3 + 
    mo_hadtotscore_grt3_imp + 
    mo_tob_gr_anyt_yn_n2 +
    ch_tabacco_passive_up_to_Y1 + 
    ch_care_main_12m_opt2_2c,
  data = bdd_final_imp_1)

## ch_briefpemocontrol_y3 ----
results %>% filter(p.value<0.05) %>% filter(outcome == "ch_briefpemocontrol_y3") %>% select(term)

model_reduit_briefpemocontrol <- lm(
  ch_briefpemocontrol_y3 ~
    ch_age_SRS_BRIEFP_Y3 + 
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
    home_total_y3 + 
    mo_hadtotscore_grt3_imp + 
    mo_tob_gr_anyt_yn_n2 +
    ch_tabacco_passive_up_to_Y1 + 
    ch_care_main_12m_opt2_2c,
  data = bdd_final_imp_1)

model_complet_briefpemocontrol <- lm(
  ch_briefpemocontrol_y3 ~
    ch_feces_rel_p4_std_Y1_10 +    
    ch_feces_rel_g2_imp_log_std_Y1 + 
    ch_feces_rel_g8_imp_log_std_Y1 +
    ch_feces_rel_g36_imp_log_std_Y1 +
    ch_feces_rel_g43_imp_log_std_Y1 +
    ch_feces_rel_g45_imp_log_std_Y1 +
    ch_feces_rel_g47_imp_log_std_Y1 +
    ch_age_SRS_BRIEFP_Y3 + 
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
    home_total_y3 + 
    mo_hadtotscore_grt3_imp + 
    mo_tob_gr_anyt_yn_n2 +
    ch_tabacco_passive_up_to_Y1 + 
    ch_care_main_12m_opt2_2c,
  data = bdd_final_imp_1)


## ch_briefpworkmemo_y3 ----
results %>% filter(p.value<0.05) %>% filter(outcome == "ch_briefpworkmemo_y3") %>% select(term)

model_reduit_briefpworkmemo <- lm(
  ch_briefpworkmemo_y3 ~
    ch_age_SRS_BRIEFP_Y3 + 
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
    home_total_y3 + 
    mo_hadtotscore_grt3_imp + 
    mo_tob_gr_anyt_yn_n2 +
    ch_tabacco_passive_up_to_Y1 + 
    ch_care_main_12m_opt2_2c,
  data = bdd_final_imp_1)

model_complet_briefpworkmemo <- lm(
  ch_briefpworkmemo_y3 ~
    ch_feces_rel_p3_std_Y1_10 +     
    ch_feces_rel_g1_imp_log_std_Y1 + 
    ch_feces_rel_g47_imp_log_std_Y1 +
    ch_feces_rel_g90_imp_log_std_Y1 +
    ch_age_SRS_BRIEFP_Y3 + 
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
    home_total_y3 + 
    mo_hadtotscore_grt3_imp + 
    mo_tob_gr_anyt_yn_n2 +
    ch_tabacco_passive_up_to_Y1 + 
    ch_care_main_12m_opt2_2c,
  data = bdd_final_imp_1)

## ch_briefpplan_y3 ----
results %>% filter(p.value<0.05) %>% filter(outcome == "ch_briefpplan_y3") %>% select(term)

model_reduit_briefpplan <- lm(
  ch_briefpplan_y3 ~
    ch_age_SRS_BRIEFP_Y3 + 
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
    home_total_y3 + 
    mo_hadtotscore_grt3_imp + 
    mo_tob_gr_anyt_yn_n2 +
    ch_tabacco_passive_up_to_Y1 + 
    ch_care_main_12m_opt2_2c,
  data = bdd_final_imp_1)

model_complet_briefpplan <- lm(
  ch_briefpplan_y3 ~
    ch_feces_rel_p1_std_Y1_10 +      
    ch_feces_rel_p3_std_Y1_10 +     
    ch_feces_rel_g29_imp_log_std_Y1 +
    ch_age_SRS_BRIEFP_Y3 + 
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
    home_total_y3 + 
    mo_hadtotscore_grt3_imp + 
    mo_tob_gr_anyt_yn_n2 +
    ch_tabacco_passive_up_to_Y1 + 
    ch_care_main_12m_opt2_2c,
  data = bdd_final_imp_1)

## ch_WPPSI_verbal_comprehension_cor_Y3 ----
results %>% filter(p.value<0.05) %>% filter(outcome == "ch_WPPSI_verbal_comprehension_cor_Y3") %>% select(term)

model_reduit_WPPSI_verbal_comprehension <- lm(
  ch_WPPSI_verbal_comprehension_cor_Y3 ~
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
    home_total_y3 + 
    mo_hadtotscore_grt3_imp + 
    mo_tob_gr_anyt_yn_n2 +
    ch_tabacco_passive_up_to_Y1 + 
    ch_care_main_12m_opt2_2c,
  data = bdd_final_imp_1)

model_complet_WPPSI_verbal_comprehension <- lm(
  ch_WPPSI_verbal_comprehension_cor_Y3 ~
    ch_feces_rel_g11_imp_log_std_Y1 +
    ch_feces_rel_g20_imp_log_std_Y1 +
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
    home_total_y3 + 
    mo_hadtotscore_grt3_imp + 
    mo_tob_gr_anyt_yn_n2 +
    ch_tabacco_passive_up_to_Y1 + 
    ch_care_main_12m_opt2_2c,
  data = bdd_final_imp_1)

## ch_WPPSI_visuospatiale_cor_Y3 ----
results %>% filter(p.value<0.05) %>% filter(outcome == "ch_WPPSI_visuospatiale_cor_Y3") %>% select(term)

## ch_WPPSI_work_memory_cor_Y3 ----
results %>% filter(p.value<0.05) %>% filter(outcome == "ch_WPPSI_work_memory_cor_Y3") %>% select(term)

model_reduit_WPPSI_work_memory <- lm(
  ch_WPPSI_work_memory_cor_Y3 ~
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
    home_total_y3 + 
    mo_hadtotscore_grt3_imp + 
    mo_tob_gr_anyt_yn_n2 +
    ch_tabacco_passive_up_to_Y1 + 
    ch_care_main_12m_opt2_2c,
  data = bdd_final_imp_1)

model_complet_WPPSI_work_memory <- lm(
  ch_WPPSI_work_memory_cor_Y3 ~
    ch_feces_rel_g27_imp_log_std_Y1 +
    ch_feces_rel_g45_imp_log_std_Y1 +
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
    home_total_y3 + 
    mo_hadtotscore_grt3_imp + 
    mo_tob_gr_anyt_yn_n2 +
    ch_tabacco_passive_up_to_Y1 + 
    ch_care_main_12m_opt2_2c,
  data = bdd_final_imp_1)

## ch_WPPSI_total_cor_Y3 ----
results %>% filter(p.value<0.05) %>% filter(outcome == "ch_WPPSI_total_cor_Y3") %>% select(term)

model_reduit_WPPSI_total <- lm(
  ch_WPPSI_total_cor_Y3 ~
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
    home_total_y3 + 
    mo_hadtotscore_grt3_imp + 
    mo_tob_gr_anyt_yn_n2 +
    ch_tabacco_passive_up_to_Y1 + 
    ch_care_main_12m_opt2_2c,
  data = bdd_final_imp_1)

model_complet_WPPSI_total <- lm(
  ch_WPPSI_total_cor_Y3 ~
    ch_feces_rel_g20_imp_log_std_Y1 +
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
    home_total_y3 + 
    mo_hadtotscore_grt3_imp + 
    mo_tob_gr_anyt_yn_n2 +
    ch_tabacco_passive_up_to_Y1 + 
    ch_care_main_12m_opt2_2c,
  data = bdd_final_imp_1)

# Assemblage ----
R2_complet_cbcl_int <- summary(model_complet_cbcl_int)$r.squared
R2_reduit_cbcl_int <- summary(model_reduit_cbcl_int)$r.squared

R2_complet_cbcl_ext <- summary(model_complet_cbcl_ext)$r.squared
R2_reduit_cbcl_ext <- summary(model_reduit_cbcl_ext)$r.squared

R2_complet_srs <- summary(model_complet_srs)$r.squared
R2_reduit_srs <- summary(model_reduit_srs)$r.squared

R2_complet_briefpinhibit <- summary(model_complet_briefpinhibit)$r.squared
R2_reduit_briefpinhibit <- summary(model_reduit_briefpinhibit)$r.squared

R2_complet_briefpshift <- summary(model_complet_briefpshift)$r.squared
R2_reduit_briefpshift <- summary(model_reduit_briefpshift)$r.squared

R2_complet_briefpemocontrol <- summary(model_complet_briefpemocontrol)$r.squared
R2_reduit_briefpemocontrol <- summary(model_reduit_briefpemocontrol)$r.squared

R2_complet_briefpworkmemo <- summary(model_complet_briefpworkmemo)$r.squared
R2_reduit_briefpworkmemo <- summary(model_reduit_briefpworkmemo)$r.squared

R2_complet_briefpplan <- summary(model_complet_briefpplan)$r.squared
R2_reduit_briefpplan <- summary(model_reduit_briefpplan)$r.squared

R2_complet_WPPSI_verbal_comprehension <- summary(model_complet_WPPSI_verbal_comprehension)$r.squared
R2_reduit_WPPSI_verbal_comprehension <- summary(model_reduit_WPPSI_verbal_comprehension)$r.squared

R2_complet_WPPSI_work_memory <- summary(model_complet_WPPSI_work_memory)$r.squared
R2_reduit_WPPSI_work_memory <- summary(model_reduit_WPPSI_work_memory)$r.squared

R2_complet_WPPSI_total <- summary(model_complet_WPPSI_total)$r.squared
R2_reduit_WPPSI_total <- summary(model_reduit_WPPSI_total)$r.squared

table_f2_Cohen <-
  list(
    ch_cbclintscore_y2 = (R2_complet_cbcl_int - R2_reduit_cbcl_int) / (1 - R2_complet_cbcl_int),
    ch_cbclextscore_y2 = (R2_complet_cbcl_ext - R2_reduit_cbcl_ext) / (1 - R2_complet_cbcl_ext),
    ch_SRStotal_y3 = (R2_complet_srs - R2_reduit_srs) / (1 - R2_complet_srs),
    ch_briefpinhibit_y3 = (R2_complet_briefpinhibit - R2_reduit_briefpinhibit) / (1 - R2_complet_briefpinhibit),
    ch_briefpshift_y3 = (R2_complet_briefpshift - R2_reduit_briefpshift) / (1 - R2_complet_briefpshift),
    ch_briefpemocontrol_y3 = (R2_complet_briefpemocontrol - R2_reduit_briefpemocontrol) / (1 - R2_complet_briefpemocontrol),
    ch_briefpworkmemo_y3 = (R2_complet_briefpworkmemo - R2_reduit_briefpworkmemo) / (1 - R2_complet_briefpworkmemo),
    ch_briefpplan_y3 = (R2_complet_briefpplan - R2_reduit_briefpplan) / (1 - R2_complet_briefpplan),
    ch_WPPSI_verbal_comprehension_cor_Y3 = (R2_complet_WPPSI_verbal_comprehension - R2_reduit_WPPSI_verbal_comprehension) / (1 - R2_complet_WPPSI_verbal_comprehension),
    ch_WPPSI_work_memory_cor_Y3 = (R2_complet_WPPSI_work_memory - R2_reduit_WPPSI_work_memory) / (1 - R2_complet_WPPSI_work_memory),
    ch_WPPSI_total_cor_Y3 = (R2_complet_WPPSI_total - R2_reduit_WPPSI_total) / (1 - R2_complet_WPPSI_total)) %>%
  unlist() %>%
  as.data.frame() %>%
  rownames_to_column(var = "outcome") %>%
  rename("f2_Cohen" = ".") %>%
  mutate(f2_Cohen = format(round(f2_Cohen, 2), digits = 2))

table_5 <- list(tidy(model_complet_cbcl_int, conf.int = TRUE) %>% mutate(outcome = "ch_cbclintscore_y2"), 
             tidy(model_complet_cbcl_ext, conf.int = TRUE) %>% mutate(outcome = "ch_cbclextscore_y2"), 
             tidy(model_complet_srs, conf.int = TRUE) %>% mutate(outcome = "ch_SRStotal_y3"), 
             tidy(model_complet_briefpinhibit, conf.int = TRUE) %>% mutate(outcome = "ch_briefpinhibit_y3"),  
             tidy(model_complet_briefpshift, conf.int = TRUE) %>% mutate(outcome = "ch_briefpshift_y3"), 
             tidy(model_complet_briefpemocontrol, conf.int = TRUE) %>% mutate(outcome = "ch_briefpemocontrol_y3"), 
             tidy(model_complet_briefpworkmemo, conf.int = TRUE) %>% mutate(outcome = "ch_briefpworkmemo_y3"), 
             tidy(model_complet_briefpplan, conf.int = TRUE) %>% mutate(outcome = "ch_briefpplan_y3"), 
             tidy(model_complet_WPPSI_verbal_comprehension, conf.int = TRUE) %>% mutate(outcome = "ch_WPPSI_verbal_comprehension_cor_Y3"), 
             tidy(model_complet_WPPSI_work_memory, conf.int = TRUE) %>% mutate(outcome = "ch_WPPSI_work_memory_cor_Y3"), 
             tidy(model_complet_WPPSI_total, conf.int = TRUE) %>% mutate(outcome = "ch_WPPSI_total_cor_Y3")) %>%
  bind_rows() %>%
  filter(term %in% c("ch_feces_rel_g19_imp_log_std_Y1", 
                     "ch_feces_rel_g21_imp_log_std_Y1", 
                     "ch_feces_rel_g36_imp_log_std_Y1", 
                     "ch_feces_rel_g45_imp_log_std_Y1", 
                     
                     "ch_feces_rel_p4_std_Y1_10",      
                     "ch_feces_rel_g1_imp_log_std_Y1", 
                     "ch_feces_rel_g13_imp_log_std_Y1",
                     "ch_feces_rel_g15_imp_log_std_Y1",
                     "ch_feces_rel_g19_imp_log_std_Y1",
                     "ch_feces_rel_g34_imp_log_std_Y1",
                     "ch_feces_rel_g39_imp_log_std_Y1",
                     "ch_feces_rel_g43_imp_log_std_Y1",
                     
                     "ch_feces_rel_p4_std_Y1_10", 
                     
                     "ch_feces_rel_p4_std_Y1_10", 
                     
                     "ch_feces_rel_g47_imp_log_std_Y1", 
                     
                     "ch_feces_rel_p4_std_Y1_10",   
                     "ch_feces_rel_g2_imp_log_std_Y1",
                     "ch_feces_rel_g8_imp_log_std_Y1",
                     "ch_feces_rel_g36_imp_log_std_Y1", 
                     "ch_feces_rel_g43_imp_log_std_Y1", 
                     "ch_feces_rel_g45_imp_log_std_Y1",
                     "ch_feces_rel_g47_imp_log_std_Y1",
                     
                     "ch_feces_rel_p3_std_Y1_10",     
                     "ch_feces_rel_g1_imp_log_std_Y1", 
                     "ch_feces_rel_g47_imp_log_std_Y1", 
                     "ch_feces_rel_g90_imp_log_std_Y1", 
                     
                     "ch_feces_rel_p1_std_Y1_10", 
                     "ch_feces_rel_p3_std_Y1_10",      
                     "ch_feces_rel_g29_imp_log_std_Y1", 
                     
                     "ch_feces_rel_g11_imp_log_std_Y1", 
                     "ch_feces_rel_g20_imp_log_std_Y1", 
                     
                     "ch_feces_rel_g27_imp_log_std_Y1", 
                     "ch_feces_rel_g45_imp_log_std_Y1", 
                     
                     "ch_feces_rel_g20_imp_log_std_Y1")) %>%
  mutate(
    estimate = format(round(estimate, 2), digits = 2),
    conf.low = format(round(conf.low, 1), digits = 1),
    conf.high = format(round(conf.high, 1), digits = 1),
    "95% CI" = paste(conf.low, conf.high, sep = ","), 
    p.value = case_when(p.value < 0.001 ~ format(round(p.value, 4), nsmall = 4),
                        p.value < 0.01 ~ format(round(p.value, 3), nsmall = 3),
                        p.value > 0.01 ~ format(round(p.value, 2), nsmall = 2))) %>%
  select(outcome, explanatory = term, Beta = estimate, "95% CI", "p.value")

table_5 <- left_join(table_5, table_f2_Cohen, by = "outcome")

corres <- data.frame(Nom = genera_var_labels, explanatory = genera_var_names, stringsAsFactors = FALSE)

table_5 <- left_join(table_5, corres, by = "explanatory")
table_5 <- table_5 %>%
  mutate(
    Nom = case_when(explanatory == "ch_feces_rel_p1_std_Y1_10" ~ "Firmicutes", 
                    explanatory == "ch_feces_rel_p3_std_Y1_10" ~ "Bacteroidetes", 
                    explanatory == "ch_feces_rel_p4_std_Y1_10" ~ "Proteobacteria", 
                    .default = Nom)) %>%
  select(outcome, explanatory, Nom, Beta, "95% CI", "p.value", f2_Cohen)

write.xlsx(table_5, file = "4_output/table_5.xlsx")


rm(model_complet_cbcl_int, model_reduit_cbcl_int, R2_complet_cbcl_int, R2_reduit_cbcl_int,
   model_complet_cbcl_ext, model_reduit_cbcl_ext, R2_complet_cbcl_ext, R2_reduit_cbcl_ext,
   model_complet_srs, model_reduit_srs, R2_complet_srs, R2_reduit_srs,
   model_complet_briefpinhibit, model_reduit_briefpinhibit, R2_complet_briefpinhibit, R2_reduit_briefpinhibit,
   model_complet_briefpshift, model_reduit_briefpshift, R2_complet_briefpshift, R2_reduit_briefpshift,
   model_complet_briefpemocontrol, model_reduit_briefpemocontrol, R2_complet_briefpemocontrol, R2_reduit_briefpemocontrol,
   model_complet_briefpworkmemo, model_reduit_briefpworkmemo, R2_complet_briefpworkmemo, R2_reduit_briefpworkmemo,
   model_complet_briefpplan, model_reduit_briefpplan, R2_complet_briefpplan, R2_reduit_briefpplan,
   model_complet_WPPSI_verbal_comprehension, model_reduit_WPPSI_verbal_comprehension, R2_complet_WPPSI_verbal_comprehension, R2_reduit_WPPSI_verbal_comprehension,
   model_complet_WPPSI_work_memory, model_reduit_WPPSI_work_memory, R2_complet_WPPSI_work_memory, R2_reduit_WPPSI_work_memory,
   model_complet_WPPSI_total, model_reduit_WPPSI_total, R2_complet_WPPSI_total, R2_reduit_WPPSI_total)

