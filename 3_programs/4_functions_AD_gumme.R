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

## Heatmap de correlation 
heatmap_cor <- function(cormat) {                  # cormat = df avec seulement les variables à mettre dans l'heatmap + des noms raccourcis
  
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  
  reorder_cormat <- function(cormat){              # Utiliser la corrélation entre les variables comme mesure de distance
    dd <- as.dist((1-cormat)/2)
    hc <- hclust(dd)
    cormat <-cormat[hc$order, hc$order]
  }
  
  cormat <- round(cor(cormat, 
                      use = "pairwise.complete.obs", 
                      method = "spearman"), 2)
  cormat <- reorder_cormat(cormat)                 # réordonner les coef de cor
  upper_tri <- get_upper_tri(cormat)               # obtenir que le triangle sup
  melted_cormat <- reshape2::melt(upper_tri, na.rm = TRUE)   # passer en df long rapidement 
  
  heatmap <-                                       # heatmap
    ggplot(melted_cormat, aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(
      low = "blue",
      high = "red",
      mid = "white",
      midpoint = 0,
      limit = c(-1, 1),
      space = "Lab",
      name = "Spearman\nCorrelation"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      size = 12,
      hjust = 1
    )) +
    coord_fixed() +
    geom_text(aes(Var2, Var1, label = value),
              color = "black",
              size = 4) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(1, 0),
      legend.position = c(0.4, 0.7),
      legend.direction = "horizontal"
    ) +
    guides(fill = guide_colorbar(
      barwidth = 7,
      barheight = 1,
      title.position = "top",
      title.hjust = 0.5
    ))
  
  return(heatmap)
}


## Visualiser les outliers d'une variable + son ident 
outliers <- function(data, var) {
  data %>% 
    select(ident, {{var}}) %>% 
    filter({{var}} %in% boxplot.stats({{var}})$out) %>% 
    arrange({{var}}) %>% 
    kable()
}


## Visualiser des corrélations entre deux groupes de variables identiques (ex: polluant i_cor VS polluants i_cor_sg)
table_cor_sg <- function(data, exposure_vec, exposure_vec_sg) {
  
  bdd_cor_sg <- data %>% 
    select(all_of({{exposure_vec}}), 
           all_of({{exposure_vec_sg}}))
  
  table_cor_sg <- 
    round(cor(bdd_cor_sg, 
              use = "pairwise.complete.obs", 
              method = "spearman"), 
          2) %>%
    as.data.frame() %>% 
    select(all_of({{exposure_vec}})) %>%
    t() %>%
    as.data.frame() %>%
    select(all_of({{exposure_vec_sg}})) %>%
    as.matrix()
  
  colnames(table_cor_sg) <- colnames(table_cor_sg) %>%
    str_replace_all(  
      c("mo_" = "",
        "ch_" = "",
        "_total_i_cor_sg_" = " ",
        "_i_cor_sg_" = " "))
  rownames(table_cor_sg) <- rownames(table_cor_sg) %>%
    str_replace_all(
      c("mo_" = "",
        "ch_" = "",
        "_total_i_cor_" = " ",
        "_i_cor_" = " "))
  table_cor_sg <- table_cor_sg %>%
    melt(na.rm = TRUE) %>%  # passer en df long rapidement 
    filter(Var1 == Var2)
  
  return(table_cor_sg)
}


## Visualiser des corrélations entre une variable et un groupe de variables (ex: xx_pool_sg_xx VS polluants)
table_cor <- function(data, var, vars){
  bdd_cor <- data %>% 
    select(all_of({{vars}}), 
           {{var}})
  table_cor <- 
    round(cor(bdd_cor, 
              method = "spearman",
              use = "pairwise.complete.obs"), 2) %>%
    as.data.frame() %>% 
    select(-{{var}}) %>%
    t()%>%
    as.data.frame()%>%
    select({{var}}) 
  return(table_cor)
}



heatmap_cor_pairwise <- function(data, vars_1, vars_2){
  
  
  bdd_cormat <- data %>% select(all_of({{vars_1}}), all_of({{vars_2}}))
  
  cormat <- round(cor(bdd_cormat, 
                      use = "pairwise.complete.obs", 
                      method = "spearman"), 2)
  
  cormat <- cormat %>% 
    as.data.frame() %>% 
    select(all_of({{vars_1}})) %>% 
    t() %>%
    as.data.frame() %>%
    select(all_of({{vars_2}})) %>%
    as.matrix()
  
  colnames(cormat) <- colnames(cormat) %>%
    str_replace_all(  
      c("mo_" = "",
        "ch_" = "",
        "_total_i_cor_sg_" = " ",
        "_i_cor_sg_" = " "))
  rownames(cormat) <- rownames(cormat) %>%
    str_replace_all(
      c("mo_" = "",
        "ch_" = "",
        "_total_i_cor_" = " ",
        "_i_cor_" = " "))
  
  cormat <- cormat %>%
    reorder_cormat()  %>%    # réordonner les coef de cor
    get_upper_tri() %>%      # obtenir que le triangle sup
    as.matrix()
  
  cormat_long <- 
    melt(cormat, na.rm = TRUE)  # passer en df long rapidement 
  
  
  heatmap <-                                       # faire la heatmap
    ggplot(cormat_long, aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(
      low = "blue",
      high = "red",
      mid = "white",
      midpoint = 0,
      limit = c(-1, 1),
      space = "Lab",
      name = "Spearman\nCorrelation"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      size = 12,
      hjust = 1
    )) +
    coord_fixed() +
    geom_text(aes(Var2, Var1, label = value),
              color = "black",
              size = 3) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(1, 0),
      legend.position = c(0.4, 0.7),
      legend.direction = "horizontal"
    ) +
    guides(fill = guide_colorbar(
      barwidth = 7,
      barheight = 1,
      title.position = "top",
      title.hjust = 0.5
    ))
  return(heatmap)
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

model_Y3_1 <- function(var_outcome, 
                       #var_microbiote, 
                       var_age) {
  model <- with(data = bdd_final_imp,
                exp = lm(var_outcome ~
                           
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
                           ch_care_main_12m_opt2_2c))
  return(model)
}

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

