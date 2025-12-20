# Investigating Scale Dependence in Large Marine Ecosystems

# Author(s): Reilly O'Connor & Marie Gutgesell
# Version: 2025-06-13

# Load Pkgs
library(tidyverse)
library(readxl)
library(ggplot2)
library(cowplot)
library(DirichletReg)
library(vegan)
library(janitor)
library(spaa)

source("../Scale_Dependence_Structural_Stability_LMEs/src/0 - Functions.R")

# Load Data
df_lme_sp <- read.csv("../Scale_Dependence_Structural_Stability_LMEs/data/LME Data/LMEs_TLassignmentbyspecies_LME.csv")
df_lme_covar <- read_xlsx("../Scale_Dependence_Structural_Stability_LMEs/data/LME Data/LME_data_14Oct2023.xlsx")

#Set Trophic Groupings
df_sp_tis <- df_lme_sp %>%
  filter(! Trophic.level == ".") %>%
  #filter(! Trophic.level < 3.0) %>%
  mutate(Trophic_Interval = case_when(
           Trophic.level < 3.0 ~ "2.0-2.9",
           Trophic.level >= 3.0  & Trophic.level < 3.5  ~ "3.0-3.4",
           Trophic.level >= 3.5  & Trophic.level < 4.0  ~ "3.5-3.9",
           Trophic.level >= 4.0  & Trophic.level < 4.5  ~ "4.0-4.4",
           Trophic.level >= 4.5 ~ "4.5-4.9",
           TRUE ~ NA_character_)) 

#Total number of species in LMEs Considered in final analysis
df_sp_num <- df_sp_tis %>%
  filter(! (LME == "Indonesian Sea" | LME == "South China Sea")) %>%
  reframe(unique_spp = unique(Species)) %>%
  reframe(total_spp = n())

df_lme <- df_sp_tis %>%
  group_by(LME) %>%
  mutate(total_sr = n()) %>%
  group_by(LME, Trophic_Interval) %>%
  reframe(sr = n(),
          total_sr = mean(total_sr),
          tips = sr/total_sr) %>%
  filter(total_sr < 2000) %>%
  arrange(total_sr, Trophic_Interval, LME)

unique(df_lme$LME)

df_lme_dir <- df_lme %>%
  dplyr::select(-sr) %>%
  # filter(! Trophic_Interval == "2.0-2.9") %>%
  pivot_wider(values_from = tips, names_from = Trophic_Interval, values_fill = 0) %>%
  mutate(scale_sr = scale(total_sr)[,1])

#Convert data for use with Dirichlet reg
df_lme_dir$Y <- DR_data(df_lme_dir[, c("2.0-2.9", "3.0-3.4", "3.5-3.9", "4.0-4.4", "4.5-4.9")])  # species proportions

#Fit model
mod1 <- DirichReg(Y ~ total_sr, data = df_lme_dir, control = list(iterlim = 6000, tol1 = 1e-6, tol2 = 1e-12), verbosity = 1)
summary(mod1)

#Generate predictions from data
newdata <- data.frame(total_sr = seq(min(df_lme_dir$total_sr), max(df_lme_dir$total_sr), length.out = 100))
pred <- predict(mod1, newdata = newdata, type = "response")

df_pred <- as.data.frame(pred)
df_pred$total_sr <- newdata$total_sr
df_long <- pivot_longer(df_pred, cols = -total_sr, names_to = "TIP", values_to = "Proportion")

df_long$TIP <- factor(df_long$TIP, levels=c('V1', 'V2', 'V3', 'V4', 'V5'))
df_lme$Trophic_Interval <- factor(df_lme$Trophic_Interval, levels=c('4.5-4.9', '4.0-4.4', '3.5-3.9', '3.0-3.4', '2.0-2.9'))


gg_lme_tips <- ggplot(df_long, aes(x = total_sr, y = Proportion, color = TIP)) +
  geom_line(linewidth = 2) +
  geom_point(data = df_lme, aes(x = total_sr, y = tips, fill = Trophic_Interval, color = NA), shape = 21, stroke = 0.5, alpha = 0.5, size = 2, color = "black") +
  labs(x = "Total Species Richness", y = "Estimated Trophic Interval Proportions", color = "Trophic Interval") +
  scale_fill_manual(values = rev(c("#7f7f7f", "#4245c4", "#f78c2f",  "#649f4d", "#b61790"))) +
  scale_color_manual(values = c("#7f7f7f", "#4245c4", "#f78c2f",  "#649f4d", "#b61790")) +
  theme_bw(base_size = 16) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  theme(legend.position = "none",
        text = element_text(family = "Arial"))

gg_lme_tips


##### Percent change per ~ 20-100 species ######
#Empty list to store ze Results
diff_list <- list()

#Window min, max n length
windows <- seq(100, 1900, by = 25)

#Loop over windows
for (start_val in windows) {
  end_val <- start_val + 25
  
  #Make new data for prediction
  newdata_perc <- data.frame(total_sr = seq(start_val, end_val, length.out = 100))
  pred_perc <- predict(mod1, newdata = newdata_perc, type = "response")
  
  #Re-format predictions
  df_pred_perc <- as.data.frame(pred_perc)
  df_pred_perc$total_sr <- newdata_perc$total_sr
  df_long_perc <- tidyr::pivot_longer(df_pred_perc, cols = -total_sr, 
                                      names_to = "TIP", values_to = "Proportion")
  
  #Calculate percentage change in the species  window
  df_prop_diff <- df_long_perc %>% 
    group_by(TIP) %>%
    summarise(
      min_sr = min(total_sr),
      max_sr = max(total_sr),
      perc_prop_at_min = 100 * Proportion[which.min(total_sr)],
      perc_prop_at_max = 100 * Proportion[which.max(total_sr)],
      perc_diff = perc_prop_at_max - perc_prop_at_min,
      .groups = "drop"
    ) %>%
    mutate(window_start = start_val, window_end = end_val)
  
  # Store in list
  diff_list[[length(diff_list) + 1]] <- df_prop_diff
}

#Combine into a single dataframe
df_all_windows <- bind_rows(diff_list)
str(df_all_windows)

df_all_windows_mean <- df_all_windows %>%
  group_by(TIP) %>%
  summarise(mean_perc_diff = mean(perc_diff),
            sd_perc_diff = sd(perc_diff),
            se_perc_diff = standard_error(perc_diff)) %>%
  mutate(Trophic_Interval = case_when(
    TIP == "V1" ~ "2.0-2.9",
    TIP == "V2" ~ "3.0-3.4",
    TIP == "V3" ~ "3.5-3.9",
    TIP == "V4" ~ "4.0-4.4",
    TIP == "V5" ~ "4.5-2.9",
    TRUE ~ NA_character_))

gg_lme_tips_perc <- ggplot(df_all_windows_mean, aes(x = mean_perc_diff, y = Trophic_Interval, color = Trophic_Interval)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_errorbarh(aes(xmin = mean_perc_diff - sd_perc_diff, xmax = mean_perc_diff + sd_perc_diff), height = 0) +
  geom_point(size = 5, shape = 21, aes(fill = Trophic_Interval), stroke = 0.5, color = "black") +
  labs(x = "Average Percent Change Per 25 Species", y = "Trophic Interval", color = "Trophic Interval") +
  scale_fill_manual(values = (c("#7f7f7f", "#4245c4", "#f78c2f",  "#649f4d", "#b61790"))) +
  scale_color_manual(values = c("#7f7f7f", "#4245c4", "#f78c2f",  "#649f4d", "#b61790")) +
  theme_bw(base_size = 16) +
  scale_x_continuous(limits = c(-3, 3), breaks = scales::pretty_breaks()) +
  theme(legend.position = "none",
        text = element_text(family = "Arial"))

gg_lme_tips_perc

##### Trophic Diversity Pyramid #####
df_sp_long_mean <- df_lme %>%
  group_by(Trophic_Interval) %>%
  reframe(TIP_mean = mean(tips, na.rm = T),
          TIP_sd = sd(tips, na.rm = T),
          TIP_se = standard_error(tips))

df_pyramid_sp <- df_sp_long_mean %>%
  mutate(Side = "Right") %>%  # Assign one side as "Right"
  bind_rows(df_sp_long_mean %>% mutate(TIP_mean = -TIP_mean, Side = "Left"))  # Mirror for "Left" side

df_pyramid_sp$Trophic_Interval <- ordered(df_pyramid_sp$Trophic_Interval,
                                          levels = c("2.0-2.9", "3.0-3.4", "3.5-3.9",  "4.0-4.4", "4.5-4.9"))
str(df_pyramid_sp)

#Plot the mirrored pyramid of richness
gg_pyramid_sp <- ggplot(df_pyramid_sp, aes(x = TIP_mean, y = Trophic_Interval, fill = Trophic_Interval)) +
  geom_col(width = 0.8, color = "black") +  # Creates bars for both sides
  geom_errorbar(aes(xmin = TIP_mean - TIP_se / 2, xmax = TIP_mean + TIP_se / 2), 
                width = 0.2, color = "black") +  # Error bars, symmetric
  labs(x = "Propotion of Species Richness", y = "Trophic Interval") +
  theme_bw(base_size = 16) +
  scale_fill_manual(values = c("#7f7f7f", "#4245c4", "#f78c2f",  "#649f4d", "#b61790")) +
  theme(
    legend.position = "none", 
    text = element_text(family = "Arial"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_blank(), 
    strip.text = element_text(face = "bold")
  )

gg_pyramid_sp

gg_pyramid_sr <- plot_grid(gg_pyramid_sp, gg_lme_tips,
                           nrow = 1, align = "hv")

gg_pyramid_sr

gg_pyramid_sr_perc <- plot_grid(gg_pyramid_sp, gg_lme_tips, gg_lme_tips_perc,
                           nrow = 1, align = "hv")

gg_pyramid_sr_perc

ggsave("../Scale_Dependence_Structural_Stability_LMEs/Figures/Figure 2 - TIPs Pyradmid - LME Richness TIPs.jpeg", plot = gg_pyramid_sr_perc, dpi = 300, width = 12, height = 4)


##### Species Turnover (Bray Curtis Similarity LMEs) #####
s_df <- df_lme_sp %>%
  dplyr::select(Species, LME) %>%
  filter(! (LME == "Indonesian Sea" | LME == "South China Sea")) %>%
  mutate(
    Species  = trimws(Species),
    LME      = trimws(LME),
    presence = 1L
  ) %>%
  distinct() %>%
  pivot_wider(
    names_from  = LME,
    values_from = presence,
    values_fill = list(presence = 0), 
    values_fn   = list(presence = max) 
  )

#Reformat matrix so is le site by species
s_matrix <- t(s_df)
s_matrix <- as.data.frame(s_matrix)
s_matrix <- row_to_names(s_matrix, row_number = 1)

#Ensure values are all numeric
s_matrix[] <- lapply(s_matrix, function(x) as.numeric(as.character(x)))
str(s_matrix)

#Calculate B-C Dissimilarity
bc_spatial <- vegdist(s_matrix, method = "bray")
bc_spatial_matrix <- as.matrix(bc_spatial)

#transform to Similarity
bc_spatial_similarity <- (1 - bc_spatial_matrix)*100

#Calculate # of LMEs in each 10% BC similarity

#Reorganize matrix
bc_spatial_similarity1<- as.dist(bc_spatial_similarity)
bc_spatial_similarity_list <- dist2list(bc_spatial_similarity1)%>%
  unite(Pair, col, row) %>%
  rename(bc_sim = "value")

#Filter out only 1 of each pair and remove same site pairs
bc_spatial_similarity_df <- bc_spatial_similarity_list %>%
  separate(Pair, into = c("Site1", "Site2"), sep = "_") %>%
  mutate(
    Site1 = pmin(Site1, Site2),
    Site2 = pmax(Site1, Site2)
  ) %>%
  filter(Site1 != Site2) %>%
  mutate(Pair = paste(Site1, Site2, sep = "_")) %>%
  distinct(Pair, .keep_all = TRUE) %>%
  dplyr::select(Pair, bc_sim) 

#Add column for each 10% bin category
bc_spatial_similarity_df <- bc_spatial_similarity_df %>%
  mutate(category = case_when(
    bc_sim >=79 ~ 79,
    bc_sim >=69 & bc_sim <79 ~ 69,
    bc_sim >=59 & bc_sim <69 ~ 59,
    bc_sim >=49 & bc_sim <59 ~ 49,
    bc_sim >=39 & bc_sim <49 ~ 39,
    bc_sim >=29 & bc_sim <39 ~ 29,
    bc_sim >=19 & bc_sim <29 ~ 19,
    bc_sim >=9 & bc_sim <19 ~ 9,
    bc_sim >=0.9 & bc_sim <9 ~ 0.9,
    bc_sim >=0.01 & bc_sim <0.9 ~ 0.01,
    bc_sim >= 0.0 ~ 0,
  ))

#Number of lmes per category for supp table
bc_summary <- bc_spatial_similarity_df %>%
  group_by(category) %>%
  count()

