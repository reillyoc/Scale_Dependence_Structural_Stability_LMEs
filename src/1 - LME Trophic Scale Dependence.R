# Investigating Scale Dependence in Large Marine Ecosystems

# Author(s): Reilly O'Connor
# Version: 2025-06-13

# Load Pkgs
library(tidyverse)
library(readxl)
library(ggplot2)
library(cowplot)
library(DirichletReg)

source("../Structural_Stability_LMEs/src/Functions.R")

# Load Data
df_lme_sp <- read.csv("../Structural_Stability_LMEs/data/LMEs_TLassignmentbyspecies_LME.csv")
df_lme_covar <- read_xlsx("../Structural_Stability_LMEs/data/LME_data_14Oct2023.xlsx")

# Alternative Trophic Groupings
# mutate(Trophic_Interval = case_when(
#   Trophic.level <= 2.8  ~ "2.0-2.8",
#   Trophic.level > 2.8  & Trophic.level <= 3.2  ~ "2.8-3.2",
#   Trophic.level > 3.2  & Trophic.level <= 3.7  ~ "3.2-3.7",
#   Trophic.level > 3.7  & Trophic.level <= 4.2  ~ "3.7-4.2",
#   Trophic.level > 4.2 ~ "4.2-5.0",
#   TRUE ~ NA_character_))


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
  filter(total_sr < 2000)

df_lme_dir <- df_lme %>%
  dplyr::select(-sr) %>%
  pivot_wider(values_from = tips, names_from = Trophic_Interval, values_fill = 0) %>%
  mutate(scale_sr = scale(total_sr)[,1])

# Convert your data to Dirichlet format
df_lme_dir$Y <- DR_data(df_lme_dir[, c("2.0-2.9", "3.0-3.4", "3.5-3.9", "4.0-4.4", "4.5-4.9")])  # species proportions

# Fit the Dirichlet regression model
mod1 <- DirichReg(Y ~ total_sr, data = df_lme_dir, control = list(iterlim = 6000, tol1 = 1e-6, tol2 = 1e-12), verbosity = 1)
summary(mod1)

newdata <- data.frame(total_sr = seq(min(df_lme_dir$total_sr), max(df_lme_dir$total_sr), length.out = 100))
pred <- predict(mod1, newdata = newdata, type = "response")

df_pred <- as.data.frame(pred)
df_pred$total_sr <- newdata$total_sr
df_long <- pivot_longer(df_pred, cols = -total_sr, names_to = "TIP", values_to = "Proportion")

gg_lme_tips <- ggplot(df_long, aes(x = total_sr, y = Proportion, color = TIP)) +
  geom_line(linewidth = 2) +
  geom_point(data = df_lme, aes(x = total_sr, y = tips, fill = Trophic_Interval, color = NA), shape = 21, stroke = 0.5, alpha = 0.5, size = 2, color = "black") +
  labs(x = "Total Species Richness", y = "Estimated Trophic Interval Proportions", color = "Trophic Interval") +
  scale_fill_manual(values = rev(c("#7f7f7f", "#4245c4", "#f78c2f",  "#649f4d", "#b61790"))) +
  scale_color_manual(values = c("#7f7f7f", "#4245c4", "#f78c2f",  "#649f4d", "#b61790")) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        text = element_text(family = "Arial"))

gg_lme_tips


##### Percent change per ~ 20-100 species ######
# Create an empty list to store results
diff_list <- list()

# Define your windows
windows <- seq(100, 1900, by = 25)

# Loop over each window
for (start_val in windows) {
  end_val <- start_val + 25
  
  # Make new data for prediction
  newdata_perc <- data.frame(total_sr = seq(start_val, end_val, length.out = 100))
  pred_perc <- predict(mod1, newdata = newdata_perc, type = "response")
  
  # Reshape predictions
  df_pred_perc <- as.data.frame(pred_perc)
  df_pred_perc$total_sr <- newdata_perc$total_sr
  df_long_perc <- tidyr::pivot_longer(df_pred_perc, cols = -total_sr, 
                                      names_to = "TIP", values_to = "Proportion")
  
  # Calculate percentage point change in the window
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

# Combine into one dataframe
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
  theme_bw(base_size = 14) +
  xlim(-3, 3) +
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

# Plot the mirrored pyramid
gg_pyramid_sp <- ggplot(df_pyramid_sp, aes(x = TIP_mean, y = Trophic_Interval, fill = Trophic_Interval)) +
  geom_col(width = 0.8, color = "black") +  # Creates bars for both sides
  geom_errorbar(aes(xmin = TIP_mean - TIP_se / 2, xmax = TIP_mean + TIP_se / 2), 
                width = 0.2, color = "black") +  # Error bars, symmetric
  labs(x = "Propotion of Species Richness", y = "Trophic Interval") +
  theme_bw(base_size = 14) +
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


ggsave("../Structural_Stability_LMEs/figures/Figure 2 - TIPs Pyradmid - LME Richness TIPs.jpeg", plot = gg_pyramid_sr_perc, dpi = 300, width = 15, height = 5)


#ACF Plots for TIPs
acf_df <- df_sp_tips_sub %>%
  group_by(Region, Trophic_interval) %>%
  summarise(acf = list(acf(tips, plot = FALSE, lag.max = 10))) %>%
  mutate(acf_df = map(acf, ~data.frame(lag = .x$lag, acf = .x$acf))) %>%
  dplyr::select(-acf) %>%
  unnest(acf_df)

ci_df <- df_sp_tips_sub %>%
  group_by(Region, Trophic_interval) %>%
  summarise(n = n(),
            ci = 1.96 / sqrt(n))

ggplot(acf_df, aes(x = as.factor(lag), y = acf)) +  
  geom_hline(data = ci_df, aes(yintercept = ci), linetype = "dashed", color = "blue", alpha = 0.5) +
  geom_hline(data = ci_df, aes(yintercept = -ci), linetype = "dashed", color = "blue", alpha = 0.5) +
  geom_hline(yintercept = 0, color = "black", alpha = 0.5) +
  geom_col() +
  facet_wrap(~Region + Trophic_interval) +
  theme_classic() +
  labs(x = "Lag", y = "ACF")

