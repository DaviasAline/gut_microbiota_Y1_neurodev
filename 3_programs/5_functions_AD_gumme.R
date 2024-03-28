## Fonctions Aline
## GUMME project

# Chargement des packages ----
library(tidyverse)
library(haven)
library(reshape2)
library(GGally)
library(gtsummary)
library(summarytools)
library(patchwork)
library(ggpubr)
library(grid)
library(questionr)
library(Hmisc)
library(rmarkdown)
library(knitr)
library(labelled)
library(distill)
library(rmdformats)
library(parameters)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(grDevices)
library(lazyeval)
library(mice)
library(glue)
library(car)
library(expss)
theme_gtsummary_language("en", decimal.mark = ".", big.mark = " ")
theme_gtsummary_compact(set_theme = TRUE)
library(bookdown)
library(performance)
library(see)


# Fonctions descriptives -----
verif_distrib <- function(data, var) {                                               # data = df large dans lequel récuperer les variables                                                                                         # vars = vecteur de noms de variables (toutes numériques)
  
  boxplot <- data %>%
    ggplot() +
    aes(x = "", y = {{var}}) +
    geom_boxplot(shape = "circle", fill = "#112446") +
    coord_flip() +
    theme_bw() +
    theme(axis.title = element_blank()) 
  densityplot <- data %>%
    ggplot() +
    aes(x = {{var}}) +
    geom_density(fill = "lightgray") +
    theme_bw() +
    stat_overlay_normal_density(color = "red", linetype = "dashed") +
    theme(axis.title = element_blank())
  qqnorm <- data %>%
    ggplot(aes(sample = {{var}})) +
    stat_qq() +
    stat_qq_line()+
    theme_bw() + 
    theme(axis.title = element_blank())
  results <- boxplot + densityplot + qqnorm
  
  return(results)
}

## Tableau descriptif variables numériques
descrip_num <- function(data, vars) {                                           
  data %>%
    select(all_of({{vars}})) %>%
    tbl_summary(
      missing = "no", 
      type = list(where(is.numeric)~ "continuous"), 
      statistic = all_continuous()~ "{min}/{p25}/{median}/{mean}/{p75}/{max}/{N_nonmiss}", 
      digits = list(all_continuous() ~ c(0, 0, 0, 1, 0, 0, 0))
    ) %>%
    bold_labels() %>% 
    as_gt()  %>% 
    as.data.frame()%>% 
    select(variable, label, stat_0) %>%
    separate(
      col = stat_0, 
      into = c("Min", "Q1", "Median", "Mean", "Q3", "Max", "N"), 
      sep = "/", 
      remove = TRUE) %>%
    rename(
      #"Variable names" = variable, 
      "Variable labels" = label)
  #%>% kable()
}

## Tableau comparaison d'effectifs
comp_effectifs <- function(data, vars_col1, vars_col2, name_col1, name_col2){ 
  table_col1 <- data %>% select(all_of({{vars_col1}})) 
  table_col2 <- data %>% select(all_of({{vars_col2}}))
  colnames(table_col2) <- colnames(table_col1)
  
  table_col1 <- table_col1 %>% tbl_summary()
  table_col2 <- table_col2 %>% tbl_summary()
  comp <- tbl_merge(
    tbls = list(table_col1, table_col2), 
    tab_spanner = c({{name_col1}}, {{name_col2}}))
  return(comp)
}


## Scatterplots
scatterplot <- function(data, outcome, vars, ncol) {
  data_long <-
    data %>%
    select(ident, all_of({{vars}}))
  var_label(data_long[, vars]) <- NULL
  data_long[, vars] <- lapply(data_long[, vars], as.numeric)
  data_long <- data_long %>%
    pivot_longer(cols = -ident,
                 names_to = "Variable",
                 values_to = "Value")
  bdd_outcome <-
    data %>% 
    select(ident, {{outcome}})
  data_long <- 
    left_join(data_long,
              bdd_outcome,
              by = "ident") %>%
    rename(Outcome = {{outcome}})
  
  scatterplot <- data_long %>%
    ggplot() +
    aes(x = Outcome, y = Value) +
    geom_point(
      shape = "circle",
      size = 1.55,
      colour = "#112446") +
    labs(y ="") +
    theme_bw() +
    facet_wrap(~Variable, scales = "free", ncol = ncol)
  return(scatterplot)
  
}

## Boxplots
boxplot <- function(data, vars, ncol) {                                               # data = df large dans lequel récuperer les variables                                                                                         # vars = vecteur de noms de variables (toutes numériques)
  data_long <-                                                     
    data %>% 
    select(ident, all_of({{vars}}))
  var_label(data_long[, vars]) <- NULL
  data_long[, vars] <- lapply(data_long[, vars], as.numeric) 
  data_long <- data_long %>% 
    pivot_longer(cols = -ident, names_to = "Variable", values_to = "Value")
  
  boxplot <- data_long %>%
    ggplot() +
    aes(x = "", y = Value) +
    geom_boxplot(shape = "circle", fill = "#112446") +
    coord_flip() +
    theme_bw() +
    theme(axis.title = element_blank()) +
    facet_wrap(~Variable, scales = "free", ncol = ncol)
  
  return(boxplot)
}


## Histogrammes
histogram <- function(data, vars, ncol) {                                             # data = df large dans lequel récuperer les variables  
  data_long <-                                                                  # vars = vecteur de noms de variables (toutes numériques)
    data %>% 
    select(ident, all_of({{vars}}))
  var_label(data_long[, vars]) <- NULL
  data_long[, vars] <- lapply(data_long[, vars], as.numeric) 
  data_long <- data_long %>% 
    pivot_longer(cols = -ident, names_to = "Variable", values_to = "Value")
  
  histogram <- 
    ggplot(data_long) +
    aes(x = Value) +
    geom_histogram(bins = 30L, fill = "#112446") +
    theme_bw() +
    theme(axis.title = element_blank()) +
    facet_wrap( ~ Variable, scales = "free", ncol = ncol)
  
  return(histogram)
} 

## Barplots 
barplot <- function(data, vars) {                                               # data = df large dans lequel récuperer les variables 
  data_long <-                                                                  # vars = vecteur de noms de variables (toutes catégorielles)
    data %>% 
    select(ident, all_of({{vars}}))
  var_label(data_long[, vars]) <- NULL
  data_long[, vars] <- lapply(data_long[, vars], factor) 
  data_long <- data_long %>% 
    pivot_longer(cols = -ident, names_to = "Variable", values_to = "Value")
  
  barplot <- data_long %>%
    filter(!is.na(Value)) %>%
    ggplot() +
    aes(x = Value) +
    geom_bar(fill = "#112446") +
    coord_flip() +
    theme_bw() +
    theme(axis.title = element_blank()) +
    facet_wrap( ~ Variable, scales = "free", ncol = 3L) 
  
  return(barplot)
}

## Densityplots 
densityplot <- function(data, vars, ncol) {                                           # data = df large dans lequel récuperer les variables  
  data_long <-                                                                  # vars = vecteur de noms de variables (toutes numériques)
    data %>% 
    select(ident, all_of({{vars}}))
  var_label(data_long[, vars]) <- NULL
  data_long[, vars] <- lapply(data_long[, vars], as.numeric) 
  data_long <- data_long %>% 
    pivot_longer(cols = -ident, names_to = "Variable", values_to = "Value")
  
  densityplot <- 
    ggplot(data_long) +
    aes(x = Value) +
    geom_density(fill = "lightgray") +
    theme_bw() +
    facet_wrap( ~ Variable, scales = "free", ncol = ncol) +
    stat_overlay_normal_density(color = "red", linetype = "dashed") +
    theme(axis.title = element_blank())
  
  return(densityplot)
} 

## qqplots
qqplot <- function(data, vars, ncol) {                                           # data = df large dans lequel récuperer les variables  
  data_long <-                                                                  # vars = vecteur de noms de variables (toutes numériques)
    data %>% 
    select(ident, all_of({{vars}}))
  var_label(data_long[, vars]) <- NULL
  data_long[, vars] <- lapply(data_long[, vars], as.numeric) 
  data_long <- data_long %>% 
    pivot_longer(cols = -ident, names_to = "Variable", values_to = "Value")
  
  qqplot <- 
    ggplot(data_long, aes(sample = Value)) +
    stat_qq() +
    stat_qq_line()+
    theme_bw() +
    facet_wrap( ~ Variable, scales = "free", ncol = ncol) +
    theme(axis.title = element_blank())
  
  return(qqplot)
} 


## Visualiser les outliers d'une variable + son ident 
outliers <- function(data, var) {
  data %>% 
    select(ident, {{var}}) %>% 
    filter({{var}} %in% boxplot.stats({{var}})$out) %>% 
    arrange({{var}}) %>% 
    kable()
}



# Fonctions statistiques ----
## Construction des modèles ----
# model_Y2_1 <- function(var_outcome, 
#                        #var_microbiote, 
#                        var_age) {
#   model <- with(data = bdd_final_imp,
#                 exp = lm(var_outcome ~
#                            
#                            #var_microbiote +
#                            
#                            var_age +      # ne pas oublier de changer selon l'outcome de neuro
#                            
#                            po_w_kg_3cat + 
#                            po_he_3cat + 
#                            mo_dipl_2cat +
#                            mo_age +
#                            mo_bmi_bepr_3cat +    
#                            ch_sex +
#                            mo_par_2cat + 
#                            ch_bf_duration_till48w_4cat +
#                            po_gd +
#                            po_delmod + 
#                            ch_food_intro_Y1_3cat +
#                            mo_pets +
#                            ch_antibio_Y1_2cat +
#                            #home_total_y3 +  
#                            mo_hadtotscore_grt3_imp +   
#                            
#                            mo_tob_gr_anyt_yn_n2+
#                            ch_tabacco_passive_up_to_Y1 +
#                            ch_care_main_12m_opt2_2c))
#   return(model)
# }
# 
# model_Y2_2 <- function(var_outcome, 
#                        #var_microbiote, 
#                        var_age) {
#   model <- with(data = bdd_final_imp,
#                 exp = lm(var_outcome ~
#                            
#                            #var_microbiote +
#                            
#                            var_age +      # ne pas oublier de changer selon l'outcome de neuro
#                            
#                            po_w_kg_3cat + 
#                            po_he_3cat + 
#                            mo_dipl_2cat +
#                            mo_age +
#                            mo_bmi_bepr_3cat +    
#                            ch_sex +
#                            mo_par_2cat + 
#                            ch_bf_duration_till48w_4cat +
#                            po_gd +
#                            po_delmod + 
#                            ch_food_intro_Y1_3cat +
#                            mo_pets +
#                            ch_antibio_Y1_2cat +
#                            #home_total_y3 +  
#                            mo_hadtotscore_grt3_imp +   
#                            
#                            mo_tob_gr_anyt_yn_n2+
#                            ch_tabacco_passive_up_to_Y1 +
#                            ch_care_main_6m_12m_opt2_2c))
#   return(model)
# }

# model_Y3_1 <- function(var_outcome, 
#                        #var_microbiote, 
#                        var_age) {
#   model <- with(data = bdd_final_imp,
#                 exp = lm(var_outcome ~
#                            
#                            #var_microbiote +
#                            
#                            var_age +      # ne pas oublier de changer selon l'outcome de neuro
#                            
#                            po_w_kg_3cat + 
#                            po_he_3cat + 
#                            mo_dipl_2cat +
#                            mo_age +
#                            mo_bmi_bepr_3cat +    
#                            ch_sex +
#                            mo_par_2cat + 
#                            ch_bf_duration_till48w_4cat +
#                            po_gd +
#                            po_delmod + 
#                            ch_food_intro_Y1_3cat +
#                            mo_pets +
#                            ch_antibio_Y1_2cat +
#                            home_total_y3 +  
#                            mo_hadtotscore_grt3_imp +   
#                            
#                            mo_tob_gr_anyt_yn_n2+
#                            ch_tabacco_passive_up_to_Y1 +
#                            ch_care_main_12m_opt2_2c))
#   return(model)
# }

# model_Y3_2 <- function(var_outcome, 
#                        #var_microbiote, 
#                        var_age) {
#   model <- with(data = bdd_final_imp,
#                 exp = lm(var_outcome ~
#                            
#                            #var_microbiote +
#                            
#                            var_age +      # ne pas oublier de changer selon l'outcome de neuro
#                            
#                            po_w_kg_3cat + 
#                            po_he_3cat + 
#                            mo_dipl_2cat +
#                            mo_age +
#                            mo_bmi_bepr_3cat +    
#                            ch_sex +
#                            mo_par_2cat + 
#                            ch_bf_duration_till48w_4cat +
#                            po_gd +
#                            po_delmod + 
#                            ch_food_intro_Y1_3cat +
#                            mo_pets +
#                            ch_antibio_Y1_2cat +
#                            home_total_y3 +  
#                            mo_hadtotscore_grt3_imp +   
#                            
#                            mo_tob_gr_anyt_yn_n2+
#                            ch_tabacco_passive_up_to_Y1 +
#                            ch_care_main_6m_12m_opt2_2c))
#   return(model)
# }

# ce modèle rend les résultats des covariables sur les outcomes, univariées, multivarié 1 et 2 (sans variables microbiote)
tbl_model_covar <- function(data, var_outcome, var_age, covar_a_tester, model_1, model_2) {
  
  effectif <- data %>%                   # Création d'une colonne effectif
    filter(!is.na({{var_outcome}}))%>%
    select(all_of({{covar_a_tester}}), 
           {{var_age}}) %>%
    tbl_summary(
      missing = "no", 
      type = list(ch_antibio_Y1_2cat ~ "categorical", 
                  mo_tob_gr_anyt_yn_n2 ~ "categorical", 
                  Mo_ETS_anyT_yn1_opt ~ "categorical", 
                  ch_ETS_12m_opt36m ~ "categorical", 
                  ch_tabacco_total_Y1 ~ "categorical", 
                  ch_tabacco_passive_up_to_Y1 ~ "categorical", 
                  home_learning_y3 ~ "continuous",
                  home_language_y3 ~ "continuous",
                  home_academic_y3 ~ "continuous", 
                  home_variety_y3 ~ "continuous",
                  home_total_y3 = "continuous"), 
      statistic = all_continuous() ~ "{N_nonmiss}", 
      
      # label = list(ch_tabacco_total_Y1 = "Child perinatal exposure to tabacco up to one year", 
      #              ch_care_main_6m_opt2_2c ~ "Main mode of child care at 6 months (2 cat)", 
      #              ch_care_main_12m_opt2_2c ~ "Main mode of child care at 12 months (2 cat)", 
      #              ch_care_main_6m_opt2_3c ~ "Main mode of child care at 6 months (3 cat)", 
      #              ch_care_main_12m_opt2_3c ~ "Main mode of child care at 12 months (3 cat)", 
      #              ch_care_main_6m_12m_op2_2c ~ "Main mode of child care up to 12 months (2 cat)",
      #              ch_care_main_6m_12m_op2_3c ~ "Main mode of child care up to 12 months (3 cat)")
    ) %>%
    bold_labels()
  
  tbl_model_univ <- data %>%                      # Création des colonnes modèle simple 
    select(
      {{var_outcome}}, 
      {{var_age}},
      all_of({{covar_a_tester}})) %>%
    
    tbl_uvregression(
      method = lm ,
      y = {{var_outcome}},
      formula = "{y} ~ {x}", 
      hide_n = TRUE, 
      estimate_fun = scales::label_number(accuracy = .1, decimal.mark = "."),
      pvalue_fun = scales::label_pvalue(accuracy = .001, decimal.mark = "."),
      # label = list(ch_tabacco_total_Y1 = "Child perinatal exposure to tabacco up to one year", 
      #              ch_care_main_6m_opt2_2c ~ "Main mode of child care at 6 months (2 cat)", 
      #              ch_care_main_12m_opt2_2c ~ "Main mode of child care at 12 months (2 cat)", 
      #              ch_care_main_6m_opt2_3c ~ "Main mode of child care at 6 months (3 cat)", 
      #              ch_care_main_12m_opt2_3c ~ "Main mode of child care at 12 months (3 cat)", 
      #              ch_care_main_6m_12m_op2_2c ~ "Main mode of child care up to 12 months (2 cat)",
      #              ch_care_main_6m_12m_op2_3c ~ "Main mode of child care up to 12 months (3 cat)")
    ) %>%
    #add_global_p() %>%
    bold_labels() %>%
    bold_p(t = 0.1)
  
  tbl_model_1 <-  tbl_regression(
    model_1, 
    estimate_fun = scales::label_number(accuracy = .1, decimal.mark = "."),
    pvalue_fun = scales::label_pvalue(accuracy = .001, decimal.mark = "."),
    # label = list(
    #   ch_tabacco_total_Y1 = "Child perinatal exposure to tabacco up to one year", 
    #   ch_care_main_6m_opt2_2c ~ "Main mode of child care at 6 months (2 cat)")
  ) %>%
    #add_global_p() %>%
    bold_labels() %>%
    bold_p(t = 0.1) %>%
    add_glance_table(include = c(nobs, r.squared, adj.r.squared)) 
  
  tbl_model_2 <-  tbl_regression(
    model_2, 
    estimate_fun = scales::label_number(accuracy = .1, decimal.mark = "."),
    pvalue_fun = scales::label_pvalue(accuracy = .001, decimal.mark = "."), 
    # label = list(
    #   ch_tabacco_total_Y1 = "Child perinatal exposure to tabacco up to one year", 
    #   ch_care_main_6m_opt2_2c ~ "Main mode of child care at 6 months (2 cat)")
  ) %>%
    #add_global_p() %>%
    bold_labels() %>%
    bold_p(t = 0.1)%>%
    add_glance_table(include = c(nobs, r.squared, adj.r.squared)) 
  
  tbl_model <- tbl_merge(
    tbls = list(effectif, 
                tbl_model_univ, 
                tbl_model_1,
                tbl_model_2), 
    tab_spanner = c("", 
                    "**Univariate models**", 
                    "**Multivariate model test 1**",
                    "**Multivariate model test 2**"))%>%
    modify_table_body(~.x %>% arrange(row_type == "glance_statistic"))
  
  return(tbl_model)
}



# check_model_Y2 <- function(
    #     var_outcome, 
#     #var_microbiote, 
#     var_age, 
#     data) {
#   model <- lm(var_outcome ~
#                 
#                 #var_microbiote +
#                 
#                 var_age +      # ne pas oublier de changer selon l'outcome de neuro
#                 
#                 po_w_kg_3cat + 
#                 po_he_3cat + 
#                 mo_dipl_2cat +
#                 mo_age +
#                 mo_bmi_bepr_3cat +    
#                 ch_sex +
#                 mo_par_2cat + 
#                 ch_bf_duration_till48w_4cat +
#                 po_gd +
#                 po_delmod + 
#                 ch_food_intro_Y1_3cat +
#                 mo_pets +
#                 ch_antibio_Y1_2cat +
#                 #home_total_y3 +  
#                 mo_hadtotscore_grt3_imp +   
#                 
#                 mo_tob_gr_anyt_yn_n2+
#                 ch_tabacco_passive_up_to_Y1 +
#                 ch_care_main_12m_opt2_2c, 
#               data = data)
#   plot <- check_model(model, check = "all")
#   
#   return(plot)
# }

check_hypothese_model <- function(
    var_outcome, 
    #var_microbiote, 
    var_age, 
    data) {
  model <- lm(var_outcome ~
                
                #var_microbiote +
                
                var_age +      # ne pas oublier de changer selon l'outcome de neuro
                
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
                
                mo_tob_gr_anyt_yn_n2+
                ch_tabacco_passive_up_to_Y1 +
                ch_care_main_12m_opt2_2c, 
              data = data)
  plot <- performance::check_model(model, check = "all"
  )
  
  return(plot)
}

# verif_residu_model_Y2 <- function(
    #     var_outcome, 
#     #var_microbiote, 
#     var_age, 
#     data) {
#   model <- lm(var_outcome ~
#                 
#                 #var_microbiote +
#                 
#                 var_age +      # ne pas oublier de changer selon l'outcome de neuro
#                 
#                 po_w_kg_3cat + 
#                 po_he_3cat + 
#                 mo_dipl_2cat +
#                 mo_age +
#                 mo_bmi_bepr_3cat +    
#                 ch_sex +
#                 mo_par_2cat + 
#                 ch_bf_duration_till48w_4cat +
#                 po_gd +
#                 po_delmod + 
#                 ch_food_intro_Y1_3cat +
#                 mo_pets +
#                 ch_antibio_Y1_2cat +
#                 #home_total_y3 +  
#                 mo_hadtotscore_grt3_imp +   
#                 
#                 mo_tob_gr_anyt_yn_n2+
#                 ch_tabacco_passive_up_to_Y1 +
#                 ch_care_main_12m_opt2_2c, 
#               data = data)
#   plot <- check_normality(model)
#   plot <- plot(plot, type = "qq")
#   
#   return(plot)
# }

verif_residu_model <- function(
    var_outcome, 
    #var_microbiote, 
    var_age, 
    data) {
  model <- lm(var_outcome ~
                
                #var_microbiote +
                
                var_age +      # ne pas oublier de changer selon l'outcome de neuro
                
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
                
                mo_tob_gr_anyt_yn_n2+
                ch_tabacco_passive_up_to_Y1 +
                ch_care_main_12m_opt2_2c, 
              data = data)
  graph <- performance::check_normality(model)
  graph <- plot(graph, type = "qq")
  
  return(graph)
}

## Présenation des résultats ----
### Tableaux analyses principales ----
# ce modèle rend les résultats des variables microbiote sur les outcomes, multivarié (tableaux pour analyses principales)
model <- function(var_outcome, var_microbiote, var_age, microbiote_name){
  model <- with(data = bdd_final_imp, 
                exp = lm(var_outcome ~ 
                           var_microbiote +
                           
                           var_age +
                           
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
                           mo_tob_gr_anyt_yn_n2+
                           ch_tabacco_passive_up_to_Y1 +
                           ch_care_main_12m_opt2_2c))
  model <- model %>% 
    tbl_regression(include = var_microbiote, 
                   label = list(var_microbiote ~ microbiote_name), 
                   estimate_fun = scales::label_number(accuracy = .01, decimal.mark = "."),
                   pvalue_fun = scales::label_pvalue(accuracy = .001, decimal.mark = ".")) %>% 
    add_n() %>% 
    bold_labels() %>%
    bold_p(t=0.1)
  return(model)
}

### Tableaux analyses supplémentaires - effets des covariables ----
tbl_model_final <- function(data, var_outcome, var_age, covar_a_tester, model_1) {
  
  effectif <- data %>%                   # Création d'une colonne effectif
    filter(!is.na({{var_outcome}}))%>%
    select(all_of({{covar_a_tester}}), 
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
      all_of({{covar_a_tester}})) %>%
    
    tbl_uvregression(
      method = lm ,
      y = {{var_outcome}},
      formula = "{y} ~ {x}", 
      hide_n = TRUE, 
      estimate_fun = scales::label_number(accuracy = .1, decimal.mark = "."),
      pvalue_fun = scales::label_pvalue(accuracy = .001, decimal.mark = "."),
      # label = list(ch_tabacco_total_Y1 = "Child perinatal exposure to tabacco up to one year", 
      #              ch_care_main_6m_opt2_2c ~ "Main mode of child care at 6 months (2 cat)", 
      #              ch_care_main_12m_opt2_2c ~ "Main mode of child care at 12 months (2 cat)", 
      #              ch_care_main_6m_opt2_3c ~ "Main mode of child care at 6 months (3 cat)", 
      #              ch_care_main_12m_opt2_3c ~ "Main mode of child care at 12 months (3 cat)", 
      #              ch_care_main_6m_12m_op2_2c ~ "Main mode of child care up to 12 months (2 cat)",
      #              ch_care_main_6m_12m_op2_3c ~ "Main mode of child care up to 12 months (3 cat)")
    ) %>%
    #add_global_p() %>%
    bold_labels() %>%
    bold_p(t = 0.1)
  
  tbl_model_1 <-  tbl_regression(
    model_1, 
    estimate_fun = scales::label_number(accuracy = .1, decimal.mark = "."),
    pvalue_fun = scales::label_pvalue(accuracy = .001, decimal.mark = "."),
    # label = list(
    #   ch_tabacco_total_Y1 = "Child perinatal exposure to tabacco up to one year", 
    #   ch_care_main_6m_opt2_2c ~ "Main mode of child care at 6 months (2 cat)")
  ) %>%
    #add_global_p() %>%
    bold_labels() %>%
    bold_p(t = 0.1) %>%
    add_glance_table(include = c(nobs, r.squared, adj.r.squared)) 
  
  tbl_model <- tbl_merge(
    tbls = list(effectif, 
                tbl_model_univ, 
                tbl_model_1), 
    tab_spanner = c("", 
                    "**Univariate models**", 
                    "**Multivariate model**"))%>%
    modify_table_body(~.x %>% arrange(row_type == "glance_statistic"))
  
  return(tbl_model)
}


### Tableux analyses supplémentaires - analyses de sensibilité ----
model_sensi_home <- function(var_outcome, var_microbiote, var_age, microbiote_name){
  model <- with(data = bdd_final_imp, 
                exp = lm(var_outcome ~ 
                           var_microbiote +
                           
                           var_age +
                           
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
  model <- model %>% 
    tbl_regression(include = var_microbiote, 
                   label = list(var_microbiote ~ microbiote_name), 
                   estimate_fun = scales::label_number(accuracy = .01, decimal.mark = "."),
                   pvalue_fun = scales::label_pvalue(accuracy = .001, decimal.mark = ".")) %>% 
    add_n() %>% 
    bold_labels() %>%
    bold_p(t=0.1)
  return(model)
}

model_sensi_seuil <- function(var_outcome, var_microbiote, var_age, microbiote_name){
  model <- with(data = bdd_final_imp_sensi_seuil, 
                exp = lm(var_outcome ~ 
                           var_microbiote +
                           
                           var_age +
                           
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
                           mo_tob_gr_anyt_yn_n2+
                           ch_tabacco_passive_up_to_Y1 +
                           ch_care_main_12m_opt2_2c))
  model <- model %>% 
    tbl_regression(include = var_microbiote, 
                   label = list(var_microbiote ~ microbiote_name), 
                   estimate_fun = scales::label_number(accuracy = .01, decimal.mark = "."),
                   pvalue_fun = scales::label_pvalue(accuracy = .001, decimal.mark = ".")) %>% 
    add_n() %>% 
    bold_labels() %>%
    bold_p(t=0.1)
  return(model)
}



# Fonction pour voir la table de correlation avant de faire la correction pour comparaison multiple selon correlation entre les tests
cor_mixed <- function(data,method="spearman"){
  ## Purpose: 2 -  use the function mixedCor from psych to calculate a 
  #                correlation matrix on mixed type variables (continuous and categorical).
  #                Note: mixedCor consider categorical variable as ordered factors.
  ## Inputs: - data: dataframe (not mice) with only the variables to consider 
  ##                 (mixed type allowed)
  ## Output: - correlation matrix
  
  ## STEP 1 : TYPE OF VARIABLES
  continuous_var = which(sapply(data, class) == "numeric")
  names(continuous_var)=NULL
  # Categorical var
  categorical_var = which(sapply(data, class) == "factor")
  #  - 2 levels only
  binary_var = categorical_var[sapply(data[,categorical_var],nlevels)==2] 
  binary_var = binary_var[!is.na(binary_var)]
  names(binary_var)=NULL
  #  - More than 2 levels (but less than 8)
  poly_var = categorical_var[(sapply(data[,categorical_var],nlevels)>2 & sapply(data[,categorical_var],nlevels)<8)] %>% na.exclude()
  names(poly_var)=NULL
  
  ## STEP 2 : CORRELATION MATRIX USING MIXEDCOR FUNCTION (FROM PSYCH)
  # data converted in numeric (necessary)
  data[,] = lapply(data[,],as.numeric)
  # Correlation matrix
  cor = data  %>% 
    mixedCor(c=continuous_var,p=poly_var,d=binary_var,use="pairwise.complete.obs",method=method)%>% pluck('rho')
  return(cor)
}


# Fonction calcul du seuil de correction pour comparaison multiple ----
# Fonction pour voir la table de correlation avant de faire la correction pour compararaison multiple selon correlation entre les tests
cor_mixed <- function(data,method="spearman"){
  ## Purpose: 2 -  use the function mixedCor from psych to calculate a 
  #                correlation matrix on mixed type variables (continuous and categorical).
  #                Note: mixedCor consider categorical variable as ordered factors.
  ## Inputs: - data: dataframe (not mice) with only the variables to consider 
  ##                 (mixed type allowed)
  ## Output: - correlation matrix
  
  ## STEP 1 : TYPE OF VARIABLES
  continuous_var = which(sapply(data, class) == "numeric")
  names(continuous_var)=NULL
  # Categorical var
  categorical_var = which(sapply(data, class) == "factor")
  #  - 2 levels only
  binary_var = categorical_var[sapply(data[,categorical_var],nlevels)==2] 
  binary_var = binary_var[!is.na(binary_var)]
  names(binary_var)=NULL
  #  - More than 2 levels (but less than 8)
  poly_var = categorical_var[(sapply(data[,categorical_var],nlevels)>2 & sapply(data[,categorical_var],nlevels)<8)] %>% na.exclude()
  names(poly_var)=NULL
  
  ## STEP 2 : CORRELATION MATRIX USING MIXEDCOR FUNCTION (FROM PSYCH)
  # data converted in numeric (necessary)
  data[,] = lapply(data[,],as.numeric)
  # Correlation matrix
  cor = data  %>% 
    mixedCor(c=continuous_var,p=poly_var,d=binary_var,use="pairwise.complete.obs",method=method)%>% pluck('rho')
  return(cor)
}

# Fonction pour voir l'alpha corrigé pour la correction pour comparaison multiple 
alpha_corrected <- function(data, alpha=0.05) {
  # Purpose: Alpha-risk correction (FWER), inspired from https://www-ncbi-nlm-nih-gov.gate2.inist.fr/pmc/articles/PMC3325408/
  # Inputs
  # - data = dataset with all exposures to consider
  # - alpha = risk to correct
  # Output
  # - alpha corrected
  
  M <- ncol(data)
  
  ## FIRST STEP: CORRELATION MATRIX (see function "cor_mixed" in the path "general_functions.R")
  cor = cor_mixed(data) #handle mixed typed variables
  
  ## SECOND STEP : EIGENVALUES OF CORRELATION MATRIX
  lambdas <- base::eigen(cor)$values
  
  ## THIRD STEP: CORRECT ALPHA
  M0 <- M - sum( (lambdas > 1) * (lambdas - 1))
  M0 <- floor(M0) # always round down estimated value
  alpha_corrected <- alpha/M0
  
  return(alpha_corrected)
  
}

# Fonction pour voir le nombre de tests rééls après prise en compte de la corrélation entre les expositions 
# Peut aussi s'appliquer à la corrélation entre les outcomes si besoin 
# puis faire nombre d'expo x nombre d'outcome
# multiplier ce nombre à la p-value
# = cette correction pour comparaison multiple adapté de Li et al 2012 est une correction de Bonferroni qui prend en compte la corrélation entre les tests effectués
M0_corrected <- function(data, alpha=0.05) {
  # Purpose: Alpha-risk correction (FWER), inspired from https://www-ncbi-nlm-nih-gov.gate2.inist.fr/pmc/articles/PMC3325408/
  # Inputs
  # - data = dataset with all exposures to consider
  # - alpha = risk to correct
  # Output
  # - alpha corrected
  
  M <- ncol(data)
  
  ## FIRST STEP: CORRELATION MATRIX (see function "cor_mixed" in the path "general_functions.R")
  cor = cor_mixed(data, method = "pearson") #handle mixed typed variables
  
  ## SECOND STEP : EIGENVALUES OF CORRELATION MATRIX
  lambdas <- base::eigen(cor)$values
  
  ## THIRD STEP: CORRECT ALPHA
  M0 <- M - sum( (lambdas > 1) * (lambdas - 1))
  M0 <- floor(M0) # always round down estimated value
  alpha_corrected <- alpha/M0
  
  return(M0)
  
}


custom_pvalue_fun <- function(x) {
  sapply(x, function(p) {
    if (is.na(p)) {
      return(NA) # Retourner NA si p est NA
    } else if (p < 0.001) {
      # Pour p < 0.001, utiliser la notation scientifique pour afficher toutes les décimales
      return(format(p, scientific = TRUE))
    } else if (p >= 0.001 & p < 0.01) {
      # Pour 0.001 <= p < 0.01, afficher avec 3 décimales
      return(sprintf("%.3f", p))
    } else {
      # Pour p >= 0.01, afficher avec 2 décimales
      return(sprintf("%.2f", p))
    }
  })
}



# Fonctions de Matthieu ----
# Functions for the preprocessing of SEPAGES data
# 25/09/23
# M. Rolland

#' Fill-in Data Below Limit Of Detection (LOD)
#'
#' This function fills in data below LOD using the fill-in method by Helsel
#' (1990). It is designed to work with data that should have a normal
#' distribution.
#'
#' @param var_to_fill A numeric vector of values to fill.
#' @param lod A numeric vector of the limit of detection (LOD), should be of the
#'   same length as `var_to_fill`.
#' @param loq A numeric vector of the limit of quantification (LOQ), should be of the
#'   same length as `var_to_fill`.
#' @details The function flags values that are to be imputed and computes
#'   distribution parameters using the \code{NADA::cenros}. It then computes
#'   fill-in values by performing a random sample between 0 and \code{lod} for a
#'   normal distribution with the previously computed parameters. Finally, the
#'   function replaces values below \code{lod} by fill-in values.
#'
#' @return A numeric vector with the original values and where those are below
#'   LOD replaced by fill-in values below the limit of detection.
#'
#' @references Helsel, D.R. (1990). Less Than Obvious - Statistical Treatment of
#'   Data Below the Detection Limit. Environmental Science & Technology, 24(12),
#'   1766-1774.
#'
#' @export
#'
#' @importFrom NADA cenros mean sd
#' @importFrom msm rtnorm
fill_in <- function(var_to_fill, lod, loq = NULL){
  # define max lim
  max_lim <- ifelse(!is.null(loq), loq, lod)
  
  # flag values to impute
  censored <- var_to_fill < max_lim 
  
  # set NA to FALSE as index cannot be NA
  censored[is.na(censored)] <- FALSE
  
  # count number of values to be imputed
  n_impute <- sum(censored)
  
  # compute distribution parameters
  stats_ros <- NADA::cenros(var_to_fill, censored, forwardT = NULL) #"forwardT = NULL": no log transform
  dist_mean <- NADA::mean(stats_ros)
  dist_sd   <- NADA::sd(stats_ros)
  
  # if data between lod and loq is to be imputed
  if(!is.null(loq)){
    n_impute_loq <- sum(between(var_to_fill, lod, loq), na.rm = TRUE)
    censored_loq <- between(var_to_fill, lod, loq)
    censored_loq[is.na(censored_loq)] <- FALSE
    
    fill_in_values_loq <- msm::rtnorm(
      n     = n_impute_loq, 
      mean  = dist_mean, 
      sd    = dist_sd, 
      lower = lod, 
      upper = loq
    )
    
    # discount varioables between lod and loq imputation of data below LOD
    censored <- ifelse(censored_loq, FALSE, censored)
    n_impute <- n_impute - n_impute_loq
  }
  
  # fill in data below lod
  fill_in_values <- msm::rtnorm(
    n     = n_impute, 
    mean  = dist_mean, 
    sd    = dist_sd, 
    lower = rep(-Inf, n_impute), 
    upper = lod
  )
  
  # replace values below LOD by fill in values
  vec_filled_in <- var_to_fill
  vec_filled_in[censored] <- fill_in_values
  if(!is.null(loq)){vec_filled_in[censored_loq] <- fill_in_values_loq}
  
  return(vec_filled_in)
}



#' Apply Standardisation on Protocol Variables
#'
#' This function standardises the exposure data on selected protocol variables
#' (either categorical or continuous) as per the sepages pipeline guide. 
#' The function is designed to be called within a `mutate` call on grouped data,
#' where each group represents a different exposure.
#'
#' @param data A data.frame in tidy format containing the data to standardise, 
#'   i.e., one row per ID and per exposure.
#' @param var_to_std A character string representing the variable to standardise.
#'   This variable should be normally distributed.
#' @param protocol_vars A character vector representing potential protocol 
#'   variables on which to standardise.
#' @param covariates A character vector of names of model covariates.
#' @param folder A character string representing the folder where standardisation 
#'   regression outputs will be saved
#' @param group A character string representing the grouping variable for 
#'   exposure. Defaults to the current grouping.
#'
#' @details
#' The function selects the final protocol variables on which to standardise 
#' based on p < 0.2. It constructs a formula, computes model residuals, 
#' sets reference values for prediction, computes standardised values, and 
#' finally returns a vector with corrected values and NA for non-computable 
#' residuals, maintaining the length same as input data.
#'
#' @return
#' A numeric vector containing the corrected values after standardisation. 
#' The length of the returned vector is the same as the number of rows in 
#' the input data. Non-computable residuals are returned as `NA`.
#' 
#' @seealso
#' \code{\link[dplyr]{mutate}}, \code{\link[stats]{lm}}, \code{\link[stats]{predict}}, \code{\link[stats]{residuals}}
#'
#' @export

standardise <- function(data = pick(everything()), var_to_std, protocol_vars, 
                        covariates = NULL, folder, group = dplyr::cur_group()){
  # Select protocol variables for standardisation (p < 0.2)
  final_std_vars <- get_protocol_var(data, var_to_std, protocol_vars, covariates, folder, group) 
  
  if(length(final_std_vars) > 0){ # apply std if at least one protocol var is associated
    
    # Construct model formula with final protocol variables (p < 0.2)
    form <- stringr::str_c(var_to_std, 
                           "~", 
                           paste(final_std_vars, collapse = "+"), 
                           if(!is.null(covariates)){"+"},
                           paste(covariates, collapse = "+")) |> 
      as.formula()
    
    # Fit model
    lm1 <- lm(form, data = data)
    
    # Get betas
    betas <- lm1 |> 
      broom::tidy() |>
      tidycat::tidy_categorical(m = lm1) 
    
    # add all correction factor
    for(prot_var in final_std_vars){
      
      # get class (numeric or categorical)
      class_i <- class(data[[prot_var]])
      
      if("numeric" %in% class_i){
        
        # get beta
        beta_i <- betas |>
          filter(variable == prot_var) |>
          pull(estimate)
        
        # prepare correction factor
        data <- data |>
          mutate(
            !!str_c("correct_", prot_var) := beta_i * (.data[[prot_var]] - median(data[[prot_var]], na.rm = TRUE))
          )
        
      }else if("factor" %in% class_i){
        
        # get betas and prepare correction factor (==beta)
        beta_i <- betas |> 
          filter(variable == prot_var) |>
          select(estimate, level) |>
          rename(
            !!prot_var := level,
            !!str_c("correct_", prot_var) := estimate
          )
        
        # add betas to data frame
        data <- data |>
          left_join(beta_i, by = prot_var)
        
      }else{
        
        stop(str_c("Protocol variables need to be coded as continuous or factor, not: ", 
                   class_i, "(variable ", prot_var, ")"))
        
      }
    }
    
    # standardise
    data <- data |>
      mutate(
        val_std = .data[[var_to_std]] - rowSums(pick(starts_with("correct")), na.rm = TRUE)
      )
  }else{
    data$val_std <- data[[var_to_std]]
  }
  
  
  
  
  # Return the vector of standardised values
  return(data$val_std)
}

#' Get Protocol Variables for Standardisation
#'
#' This function identifies which among the potential protocol variables should
#' be used for standardising a given exposure. It exports outputs for the
#' different steps of the process.
#'
#' @param data A data frame containing the exposure and protocol variables.
#' @param var_to_std A character string with the name of the variable to be 
#'   standardised.
#' @param protocol_vars A character vector of names of potential protocol
#'   variables.
#' @param covariates A character vector of names of model covariates.
#' @param folder The directory folder where the output CSV files will be saved.
#' @param group The group/exposure under consideration.
#' @return A character vector containing the names of the final protocol
#'   variables to be used for standardisation.
#' 
get_protocol_var <- function(data, var_to_std, protocol_vars, covariates, folder, group){
  # Construct linear model formula
  model_formula <- str_c(var_to_std, 
                         "~", 
                         paste(protocol_vars, collapse = "+"), 
                         if(!is.null(covariates)){"+"},
                         paste(covariates, collapse = "+"))
  
  # Fit the full linear model
  lm_full <- lm(as.formula(model_formula), data = data)
  
  # Extract betas from the model and export to CSV
  betas <- broom::tidy(lm_full)
  
  # export
  filename <- str_c(paste(group, collapse = "_"), ".csv")
  write_csv(betas, file.path(folder, filename))
  
  # Perform ANOVA and get p-values
  aov_output <- car::Anova(lm_full)
  
  # export anova output
  filename <- str_c(paste(group, collapse = "_"), "_aov.csv") 
  write_csv(broom::tidy(aov_output), file.path(folder, filename))
  
  # Identify and return protocol variables with p < 0.2
  final_std_vars <- aov_output |>
    as.data.frame() |> 
    rownames_to_column(var = "term") |>
    rename(p = "Pr(>F)") |>
    filter(term %in% protocol_vars & p < 0.2) |>
    pull(term)
  
  return(final_std_vars)
}
