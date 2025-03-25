## Fig.1: Forestplot phyla ----
figure_1 <- table_multi %>% 
  mutate(
    q_value_shape = ifelse(`p-value`<0.5, "p-value<0.05", "p-value≥0.05")) %>%
  filter(`Gut microbiota parameters` %in% c("Firmicutes",
                                            "Actinobacteria", 
                                            "Bacteroidetes", 
                                            "Proteobacteria")) %>%
  mutate(Beta = as.numeric(Beta), 
         Outcome_rec = gsub("WPPSI", "WPPSI-IV", Outcome_rec), 
         Outcome_rec = gsub("SRS", "SRS-II", Outcome_rec), 
         Outcome_rec =
           fct_relevel(Outcome_rec,
                       "Total WPPSI-IV Y3", "Work memory WPPSI-IV Y3",
                       "Visuospatial WPPSI-IV Y3","Verbal comprehension WPPSI-IV Y3",
                       "Plan and organization BRIEF-P Y3", "Work memory BRIEF-P Y3",
                       "Emotional control BRIEF-P Y3","Shift BRIEF-P Y3",
                       "Inhibition BRIEF-P Y3", "Total SRS-II Y3",
                       "Externalizing CBCL Y2", "Internalizing CBCL Y2"),
         # Outcome_rec = fct_rev(Outcome_rec),
         `Gut microbiota parameters` =
           fct_relevel(`Gut microbiota parameters`, 
                       "Proteobacteria", "Bacteroidetes", "Actinobacteria", "Firmicutes" )) %>%
  ggplot(aes(x = `Gut microbiota parameters`, 
             y = Beta, 
             min = lower_CI, 
             ymax = upper_CI, 
             color = q_value_shape)) +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_pointrange(
    position = position_dodge(width = 0.5), 
    size = 0.4) +
  labs(x = "Neurodevelopement", y = "", color = "p-value") +
  theme_lucid() +
  coord_flip()  +
  facet_wrap(`Outcome_rec`~.,  ncol = 12) +
  theme(
    legend.position = "right",
    legend.box = "vertical", 
    legend.justification = "right")  +
  scale_color_manual(values = c("p-value<0.05" = "red", "p-value≥0.05" = "black"))

figure_1
ggsave("C:/Users/Aline/OneDrive - etu.univ-grenoble-alpes.fr/Documents/7. Présentations orales/17. Thèse/forest_plot_phyla_neuro.tiff", 
       figure_1, 
       device = "tiff",
       units = "cm",
       dpi = 300,
       height = 15, 
       width = 25)



## Fig.1: Forestplot phyla ----
figure_1 <- table_multi %>% 
  mutate(
    q_value_shape = ifelse(`p-value`<0.0002, "p-value<0.0002", "p-value≥0.0002")) %>%
  filter(`Gut microbiota parameters` %in% c("Firmicutes",
                                            "Actinobacteria", 
                                            "Bacteroidetes", 
                                            "Proteobacteria")) %>%
  mutate(Beta = as.numeric(Beta), 
         Outcome_rec = gsub("WPPSI", "WPPSI-IV", Outcome_rec), 
         Outcome_rec = gsub("SRS", "SRS-II", Outcome_rec), 
         Outcome_rec =
           fct_relevel(Outcome_rec,
                       "Total WPPSI-IV Y3", "Work memory WPPSI-IV Y3",
                       "Visuospatial WPPSI-IV Y3","Verbal comprehension WPPSI-IV Y3",
                       "Plan and organization BRIEF-P Y3", "Work memory BRIEF-P Y3",
                       "Emotional control BRIEF-P Y3","Shift BRIEF-P Y3",
                       "Inhibition BRIEF-P Y3", "Total SRS-II Y3",
                       "Externalizing CBCL Y2", "Internalizing CBCL Y2"),
         # Outcome_rec = fct_rev(Outcome_rec),
         `Gut microbiota parameters` =
           fct_relevel(`Gut microbiota parameters`, 
                       "Proteobacteria", "Bacteroidetes", "Actinobacteria", "Firmicutes" )) %>%
  ggplot(aes(x = `Gut microbiota parameters`, 
             y = Beta, 
             min = lower_CI, 
             ymax = upper_CI, 
             color = q_value_shape)) +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_pointrange(
    position = position_dodge(width = 0.5), 
    size = 0.4) +
  labs(x = "Neurodevelopement", y = "", color = "p-value") +
  theme_lucid() +
  coord_flip()  +
  facet_wrap(`Outcome_rec`~.,  ncol = 12) +
  theme(
    legend.position = "right",
    legend.box = "vertical", 
    legend.justification = "right")  +
  scale_color_manual(values = c("p-value<0.0002" = "red", "p-value≥0.0002" = "black"))

figure_1
ggsave("C:/Users/Aline/OneDrive - etu.univ-grenoble-alpes.fr/Documents/7. Présentations orales/17. Thèse/forest_plot_phyla_neuro.tiff", 
       figure_1, 
       device = "tiff",
       units = "cm",
       dpi = 300,
       height = 15, 
       width = 25)




## Fig.1: Forestplot phyla ----
figure_1 <- table_multi %>% 
  mutate(
    q_value_shape = ifelse(`p-value`<0.0002, "p-value<0.0002", "p-value≥0.0002")) %>%
  filter(`Gut microbiota parameters` %in% c("Firmicutes",
                                            "Actinobacteria", 
                                            "Bacteroidetes", 
                                            "Proteobacteria")) %>%
  mutate(Beta = as.numeric(Beta), 
         Outcome_rec = gsub("WPPSI", "WPPSI-IV", Outcome_rec), 
         Outcome_rec = gsub("SRS", "SRS-II", Outcome_rec), 
         Outcome_rec =
           fct_relevel(Outcome_rec,
                       "Total WPPSI-IV Y3", "Work memory WPPSI-IV Y3",
                       "Visuospatial WPPSI-IV Y3","Verbal comprehension WPPSI-IV Y3",
                       "Plan and organization BRIEF-P Y3", "Work memory BRIEF-P Y3",
                       "Emotional control BRIEF-P Y3","Shift BRIEF-P Y3",
                       "Inhibition BRIEF-P Y3", "Total SRS-II Y3",
                       "Externalizing CBCL Y2", "Internalizing CBCL Y2"),
         # Outcome_rec = fct_rev(Outcome_rec),
         `Gut microbiota parameters` =
           fct_relevel(`Gut microbiota parameters`, 
                       "Firmicutes", "Actinobacteria", "Bacteroidetes", "Proteobacteria")) %>%
  ggplot(aes(x = Outcome_rec, 
             y = Beta, 
             min = lower_CI, 
             ymax = upper_CI, 
             color = q_value_shape)) +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_pointrange(
    position = position_dodge(width = 0.5), 
    size = 0.4) +
  labs(x = "Neurodevelopement", y = "", color = "p-value") +
  theme_lucid() +
  coord_flip()  +
  facet_wrap(`Gut microbiota parameters`~.,  ncol = 4) +
  theme(
    legend.position = "right",
    legend.box = "vertical", 
    legend.justification = "right")  +
  scale_color_manual(values = c("p-value<0.0002" = "red", "p-value≥0.0002" = "black"))

figure_1
ggsave("C:/Users/Aline/OneDrive - etu.univ-grenoble-alpes.fr/Documents/7. Présentations orales/17. Thèse/forest_plot_phyla_neuro.tiff", 
       figure_1, 
       device = "tiff",
       units = "cm",
       dpi = 300,
       height = 15, 
       width = 25)




figure_3 <- table_multi %>% 
  mutate(Beta = as.numeric(Beta), 
         Outcome_rec = gsub("WPPSI", "WPPSI-IV", Outcome_rec), 
         Outcome_rec = gsub("SRS", "SRS-II", Outcome_rec), 
         Outcome_rec = 
           fct_relevel(Outcome_rec,
                       "Internalizing CBCL Y2", 
                       "Externalizing CBCL Y2",
                       "Total SRS-II Y3", 
                       "Inhibition BRIEF-P Y3", 
                       "Shift BRIEF-P Y3",
                       "Emotional control BRIEF-P Y3",
                       "Work memory BRIEF-P Y3",
                       "Plan and organization BRIEF-P Y3", 
                       "Verbal comprehension WPPSI-IV Y3", 
                       "Visuospatial WPPSI-IV Y3",
                       "Work memory WPPSI-IV Y3", 
                       "Total WPPSI-IV Y3"), 
         `Gut microbiota parameters` = 
           fct_relevel(`Gut microbiota parameters`, 
                       "Saccharibacteria genera incertae sedis", "Peptoniphilus",
                       "Granulicatella", "Anaerotruncus", "Lactococcus", "Terrisporobacter",
                       "Oscillibacter", "Haemophilus", "Erysipelotrichaceae incertae sedis",
                       "Butyricicoccus", "Dialister", "Subdoligranulum", "Intestinibacter",
                       "Klebsiella", "Eisenbergiella", "Hungatella", "Dorea", "Eggerthella",
                       "Romboutsia", "Clostridium IV", "Ruminococcus 2", "Flavonifractor",
                       "Alistipes", "Collinsella", "Parabacteroides", "Fusicatenibacter",
                       "Roseburia", "Enterobacter", "Cellulosibacter", "Enterococcus",
                       "Coprococcus", "Veillonella", "Clostridium sensu stricto", "Clostridium XVIII",
                       "Ruminococcus", "Gemmiger", "Anaerostipes", "Lachnospiracea incertae sedis",
                       "Clostridium XlVa", "Streptococcus", "Faecalibacterium", "Akkermansia",
                       "Escherichia and Shigella", "Blautia", "Bacteroides", "Bifidobacterium",
                       "Proteobacteria", "Bacteroidetes", "Actinobacteria", "Firmicutes",
                       "Shannon diversity", "Specific richness"), 
         Phyla_corres = fct_recode(Phyla_corres, "Candidatus Saccharibacteria" = "Candidatus_Saccharibacteria"), 
         Phyla_corres = fct_relevel(Phyla_corres,
                                    "Firmicutes", "Actinobacteria", "Bacteroidetes", "Proteobacteria",
                                    "Verrucomicrobia", "Candidatus Saccharibacteria"), 
         `Gut microbiota parameters` = 
           fct_recode(`Gut microbiota parameters`, 
                      "Saccharibacteria" = "Saccharibacteria genera incertae sedis")) %>%
  #filter(`p-value`<0.05) %>% 
  filter(`Gut microbiota parameters` %in% c("Firmicutes",
                                             "Actinobacteria",
                                             "Bacteroidetes",
                                             "Proteobacteria")) %>%
  ggplot(aes(x = `Gut microbiota parameters`, 
             y = Beta, 
             min = lower_CI, 
             ymax = upper_CI 
             # color = Phyla_corres,
             # label = ifelse(`p-value` < 0.05, as.character(`Gut microbiota parameters`), "")
             )) +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_pointrange(
    size = 0.4,
    position = position_dodge(width = 1.2, preserve = "total")) +
  # geom_text_repel(
  #   fontface = "italic", 
  #   nudge_x  = 0.3,
  #   aes(label = ifelse(`p-value` < 0.05, as.character(`Gut microbiota parameters`), "")) 
  #   #         fontface = "italic", 
  #   #         hjust = 0.5, 
  #   #         vjust = -0.5, 
  #   #         angle = 0, 
  #   #         size = 4, 
  #   #         position = position_dodge(width = 1.2, preserve = "total")
  # ) +
  coord_flip()  +
  facet_grid(.~Outcome_rec, 
             scales = "free_y", 
             space = "free_y", 
             switch = "y", 
             labeller = as_labeller(function(labels) {
               labels <- gsub("Emotional control BRIEF-P Y3", "Emotional control\nBRIEF-P Y3", labels)
               labels <- gsub("Work memory", "Work memory\n", labels)
               labels <- gsub("Plan and organization BRIEF-P Y3", "Plan and organization\nBRIEF-P Y3", labels)
               labels <- gsub("Verbal comprehension WPPSI-IV Y3", "Verbal comprehension\nWPPSI-IV Y3", labels)
               return(labels)
             })
  ) + 
  labs(x = "", y = "") +
  # scale_color_manual(values = c(
  #   "Firmicutes" = "black", 
  #   "Actinobacteria" = "#1b1b8e", # dark blue
  #   "Bacteroidetes" = "#1b8e1b", # dark green
  #   "Proteobacteria" = "darkred", # dark red
  #   "Verrucomicrobia" = "#6b1b8e", # dark purple
  #   "Candidatus Saccharibacteria" = "#8e6b1b")) +   # dark brown
  theme_lucid()  + 
  theme(axis.text.y = element_blank(), 
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 10), 
        panel.background = element_rect(color = "gray", fill = NULL), 
        legend.position = "bottom", 
        legend.title = element_text(face = "bold"))
figure_3


ggsave("4_output/figures/Fig.3 forest_plot_genera.tiff", 
       figure_3, 
       device = "tiff",
       units = "cm",
       dpi = 300,
       height = 20, 
       width = 34)


figure_3 <- table_multi %>% 
  mutate(Beta = as.numeric(Beta), 
         Outcome_rec = gsub("WPPSI", "WPPSI-IV", Outcome_rec), 
         Outcome_rec = gsub("SRS", "SRS-II", Outcome_rec), 
         Outcome_rec = 
           fct_relevel(Outcome_rec,
                       "Internalizing CBCL Y2", 
                       "Externalizing CBCL Y2",
                       "Total SRS-II Y3", 
                       "Inhibition BRIEF-P Y3", 
                       "Shift BRIEF-P Y3",
                       "Emotional control BRIEF-P Y3",
                       "Work memory BRIEF-P Y3",
                       "Plan and organization BRIEF-P Y3", 
                       "Verbal comprehension WPPSI-IV Y3", 
                       "Visuospatial WPPSI-IV Y3",
                       "Work memory WPPSI-IV Y3", 
                       "Total WPPSI-IV Y3"), 
         `Gut microbiota parameters` = 
           fct_relevel(`Gut microbiota parameters`, 
                       "Saccharibacteria genera incertae sedis", "Peptoniphilus",
                       "Granulicatella", "Anaerotruncus", "Lactococcus", "Terrisporobacter",
                       "Oscillibacter", "Haemophilus", "Erysipelotrichaceae incertae sedis",
                       "Butyricicoccus", "Dialister", "Subdoligranulum", "Intestinibacter",
                       "Klebsiella", "Eisenbergiella", "Hungatella", "Dorea", "Eggerthella",
                       "Romboutsia", "Clostridium IV", "Ruminococcus 2", "Flavonifractor",
                       "Alistipes", "Collinsella", "Parabacteroides", "Fusicatenibacter",
                       "Roseburia", "Enterobacter", "Cellulosibacter", "Enterococcus",
                       "Coprococcus", "Veillonella", "Clostridium sensu stricto", "Clostridium XVIII",
                       "Ruminococcus", "Gemmiger", "Anaerostipes", "Lachnospiracea incertae sedis",
                       "Clostridium XlVa", "Streptococcus", "Faecalibacterium", "Akkermansia",
                       "Escherichia and Shigella", "Blautia", "Bacteroides", "Bifidobacterium",
                       "Proteobacteria", "Bacteroidetes", "Actinobacteria", "Firmicutes",
                       "Shannon diversity", "Specific richness"), 
         Phyla_corres = fct_recode(Phyla_corres, "Candidatus Saccharibacteria" = "Candidatus_Saccharibacteria"), 
         Phyla_corres = fct_relevel(Phyla_corres,
                                    "Firmicutes", "Actinobacteria", "Bacteroidetes", "Proteobacteria",
                                    "Verrucomicrobia", "Candidatus Saccharibacteria"), 
         `Gut microbiota parameters` = 
           fct_recode(`Gut microbiota parameters`, 
                      "Saccharibacteria" = "Saccharibacteria genera incertae sedis"), 
         q_value_shape = ifelse(`p-value`<0.5, "p-value<0.05", "p-value≥0.05")) %>%
  #filter(`p-value`<0.05) %>% 
  filter(`Gut microbiota parameters` %in% c("Firmicutes",
                                            "Actinobacteria",
                                            "Bacteroidetes",
                                            "Proteobacteria")) %>%
  ggplot(aes(x = `Gut microbiota parameters`, 
             y = Beta, 
             min = lower_CI, 
             ymax = upper_CI, 
             color = p_value_shape)) +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_pointrange(
    size = 0.4,
    position = position_dodge(width = 1.2, preserve = "total")) +
  coord_flip()  +
  facet_grid(.~Outcome_rec, 
             scales = "free_y", 
             space = "free_y", 
             switch = "y", 
             labeller = as_labeller(function(labels) {
               labels <- gsub("Externalizing CBCL Y2", "Externalizing\nCBCL Y2", labels)
               labels <- gsub("Internalizing CBCL Y2", "Internalizing\nCBCL Y2", labels)
               labels <- gsub("Emotional control BRIEF-P Y3", "Emotional control\nBRIEF-P Y3", labels)
               labels <- gsub("Emotional control BRIEF-P Y3", "Emotional control\nBRIEF-P Y3", labels)
               labels <- gsub("Work memory", "Work memory\n", labels)
               labels <- gsub("Plan and organization BRIEF-P Y3", "Plan and organization\nBRIEF-P Y3", labels)
               labels <- gsub("Verbal comprehension WPPSI-IV Y3", "Verbal comprehension\nWPPSI-IV Y3", labels)
               labels <- gsub("Visuospatial WPPSI-IV Y3", "Visuospatial\nWPPSI-IV Y3", labels)
               return(labels)
             })
  ) + 
  labs(x = "", y = "", color = "p-value") +
  theme_lucid()  + 
  theme(panel.background = element_rect(color = "gray", fill = NULL), 
        legend.position = "bottom", 
        legend.title = element_text(face = "bold")) +
  scale_color_manual(values = c("p-value<0.05" = "red", "p-value≥0.05" = "black"))
figure_3

ggsave("C:/Users/Aline/OneDrive - etu.univ-grenoble-alpes.fr/Documents/7. Présentations orales/17. Thèse/forest_plot_phyla_neuro_0.05.tiff", 
       figure_3, 
       device = "tiff",
       units = "cm",
       dpi = 300,
       height = 15, 
       width = 45)




figure_3 <- table_multi %>% 
  mutate(Beta = as.numeric(Beta), 
         Outcome_rec = gsub("WPPSI", "WPPSI-IV", Outcome_rec), 
         Outcome_rec = gsub("SRS", "SRS-II", Outcome_rec), 
         Outcome_rec = 
           fct_relevel(Outcome_rec,
                       "Internalizing CBCL Y2", 
                       "Externalizing CBCL Y2",
                       "Total SRS-II Y3", 
                       "Inhibition BRIEF-P Y3", 
                       "Shift BRIEF-P Y3",
                       "Emotional control BRIEF-P Y3",
                       "Work memory BRIEF-P Y3",
                       "Plan and organization BRIEF-P Y3", 
                       "Verbal comprehension WPPSI-IV Y3", 
                       "Visuospatial WPPSI-IV Y3",
                       "Work memory WPPSI-IV Y3", 
                       "Total WPPSI-IV Y3"), 
         `Gut microbiota parameters` = 
           fct_relevel(`Gut microbiota parameters`, 
                       "Saccharibacteria genera incertae sedis", "Peptoniphilus",
                       "Granulicatella", "Anaerotruncus", "Lactococcus", "Terrisporobacter",
                       "Oscillibacter", "Haemophilus", "Erysipelotrichaceae incertae sedis",
                       "Butyricicoccus", "Dialister", "Subdoligranulum", "Intestinibacter",
                       "Klebsiella", "Eisenbergiella", "Hungatella", "Dorea", "Eggerthella",
                       "Romboutsia", "Clostridium IV", "Ruminococcus 2", "Flavonifractor",
                       "Alistipes", "Collinsella", "Parabacteroides", "Fusicatenibacter",
                       "Roseburia", "Enterobacter", "Cellulosibacter", "Enterococcus",
                       "Coprococcus", "Veillonella", "Clostridium sensu stricto", "Clostridium XVIII",
                       "Ruminococcus", "Gemmiger", "Anaerostipes", "Lachnospiracea incertae sedis",
                       "Clostridium XlVa", "Streptococcus", "Faecalibacterium", "Akkermansia",
                       "Escherichia and Shigella", "Blautia", "Bacteroides", "Bifidobacterium",
                       "Proteobacteria", "Bacteroidetes", "Actinobacteria", "Firmicutes",
                       "Shannon diversity", "Specific richness"), 
         Phyla_corres = fct_recode(Phyla_corres, "Candidatus Saccharibacteria" = "Candidatus_Saccharibacteria"), 
         Phyla_corres = fct_relevel(Phyla_corres,
                                    "Firmicutes", "Actinobacteria", "Bacteroidetes", "Proteobacteria",
                                    "Verrucomicrobia", "Candidatus Saccharibacteria"), 
         `Gut microbiota parameters` = 
           fct_recode(`Gut microbiota parameters`, 
                      "Saccharibacteria" = "Saccharibacteria genera incertae sedis"), 
         q_value_shape = ifelse(`p-value`<0.0002, "p-value<0.0002", "p-value≥0.0002")) %>%
  #filter(`p-value`<0.05) %>% 
  filter(`Gut microbiota parameters` %in% c("Firmicutes",
                                            "Actinobacteria",
                                            "Bacteroidetes",
                                            "Proteobacteria")) %>%
  ggplot(aes(x = `Gut microbiota parameters`, 
             y = Beta, 
             min = lower_CI, 
             ymax = upper_CI, 
             color = p_value_shape)) +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_pointrange(
    size = 0.4,
    position = position_dodge(width = 1.2, preserve = "total")) +
  coord_flip()  +
  facet_grid(.~Outcome_rec, 
             scales = "free_y", 
             space = "free_y", 
             switch = "y", 
             labeller = as_labeller(function(labels) {
               labels <- gsub("Externalizing CBCL Y2", "Externalizing\nCBCL Y2", labels)
               labels <- gsub("Internalizing CBCL Y2", "Internalizing\nCBCL Y2", labels)
               labels <- gsub("Emotional control BRIEF-P Y3", "Emotional control\nBRIEF-P Y3", labels)
               labels <- gsub("Emotional control BRIEF-P Y3", "Emotional control\nBRIEF-P Y3", labels)
               labels <- gsub("Work memory", "Work memory\n", labels)
               labels <- gsub("Plan and organization BRIEF-P Y3", "Plan and organization\nBRIEF-P Y3", labels)
               labels <- gsub("Verbal comprehension WPPSI-IV Y3", "Verbal comprehension\nWPPSI-IV Y3", labels)
               labels <- gsub("Visuospatial WPPSI-IV Y3", "Visuospatial\nWPPSI-IV Y3", labels)
               return(labels)
             })
  ) + 
  labs(x = "", y = "", color = "p-value") +
  theme_lucid()  + 
  theme(panel.background = element_rect(color = "gray", fill = NULL), 
        legend.position = "bottom", 
        legend.title = element_text(face = "bold")) +
  scale_color_manual(values = c("p-value<0.0002" = "red", "p-value≥0.0002" = "black"))
figure_3

ggsave("C:/Users/Aline/OneDrive - etu.univ-grenoble-alpes.fr/Documents/7. Présentations orales/17. Thèse/forest_plot_phyla_neuro_0.0002.tiff", 
       figure_3, 
       device = "tiff",
       units = "cm",
       dpi = 300,
       height = 15, 
       width = 45)

