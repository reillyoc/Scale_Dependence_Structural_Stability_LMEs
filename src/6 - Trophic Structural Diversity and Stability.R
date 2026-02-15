# Shannon Structural Diversity and Structural Variability

# Author(s): Reilly O'Connor
# Version: 2025-11-13

# Load Pkgs
library(tidyverse)
library(vegan)
library(reshape2)
library(zoo)
library(cowplot)

# load data
df_sp_tips <- read.csv("../Scale_Dependence_Structural_Stability_LMEs/Data/Output Data/richness_tips_jack2.csv", header = T) %>%
  filter(year < 2004)

#Regional Grounfish Biomass Collapses
df_gfb_windows <- tribble(
  ~Region, ~xintercept,
  "NFLS",  1990,
  # "NFLS",  1994,
  "SS",    1992,
  # "SS",    1998,
  "NEUS",  1990,
  # "NEUS",  1991,
)

df_sp_tips_t <- df_sp_tips %>%
  select(year, Region, Trophic_interval, jack2) %>%
  pivot_wider(id_cols = c(year, Region),
              names_from = Trophic_interval,
              values_from = jack2)

df_nfls_t <- df_sp_tips_t %>%
  filter(Region == "NFLS")

df_nfls_entropy <- df_nfls_t %>%
  rowwise() %>%
  mutate(shannon = diversity(c_across(`3.0-3.4`:`4.5-4.9`), index = "shannon")) %>%
  ungroup()

df_ss_t <- df_sp_tips_t %>%
  filter(Region == "SS")


df_ss_entropy <- df_ss_t %>%
  rowwise() %>%
  mutate(shannon = diversity(c_across(`3.0-3.4`:`4.5-4.9`), index = "shannon")) %>%
  ungroup()


df_neus_t <- df_sp_tips_t %>%
  filter(Region == "NEUS")


df_neus_entropy <- df_neus_t %>%
  rowwise() %>%
  mutate(shannon = diversity(c_across(`3.0-3.4`:`4.5-4.9`), index = "shannon")) %>%
  ungroup()

df_lme_entropy <- rbind(df_nfls_entropy, df_ss_entropy, df_neus_entropy)

df_lme_entropy$Region <- factor(df_lme_entropy$Region, levels=c('NFLS', 'SS', 'NEUS'))
df_gfb_windows$Region <- factor(df_gfb_windows$Region, levels=c('NFLS', 'SS', 'NEUS'))

gg_lme_entropy <- ggplot(df_lme_entropy, aes(x = year, y = shannon, color = Region)) +
  geom_vline(data = df_gfb_windows, aes(xintercept = xintercept), linetype = "dashed", alpha = 0.8) +
  facet_wrap(~ Region, scale = "free", ncol= 1) +
  geom_line(alpha = 0.8, size = 2) +
  # geom_point() +
  labs(x = "Year",
       y = "Trophic \nStructural Diversity") +
  scale_color_viridis_d() +
  scale_y_continuous(limits = c(0.90, 1.4), breaks = scales::pretty_breaks()) +
  scale_x_continuous(limits = c(1970, 2005), breaks = scales::pretty_breaks()) +
  theme_bw(base_size = 18) +
  theme(text = element_text(family = "Arial"),
        legend.position = 'none')

gg_lme_entropy

gg_tips_ssd <- plot_grid(gg_tips, gg_lme_entropy, ncol = 2, align = "hv")

gg_tips_ssd

# ggsave("../Scale_Dependence_Structural_Stability_LMEs/Figures/Figure 3 - Temporal TIPs + Structural SR.jpeg", plot = gg_tips_ssd, dpi = 300, width = 10, height = 10)


acf1_fun <- function(x) {
  if (all(is.na(x))) return(NA_real_)  # handle all-NA windows
  acf(x, lag.max = 1, plot = FALSE)$acf[2]
}


#Kvalseths CV
k_cv <- function(x) {
  if (length(x) < 2) return(NA)
  CV <- sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
  k_CV <- sqrt( (CV^2) / (1 + (CV^2)) )
  return(k_CV)
  
}

window_size <- 5

rolling_metrics <- df_lme_entropy %>%
  group_by(Region) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(across(
    .cols = c(matches("^[0-9]"), shannon),   # trophic bins + shannon
    .fns  = list(
      roll_mean = ~ rollapply(.x, width = window_size, FUN = mean, fill = NA, align = "right", partial = F),
      roll_sd   = ~ rollapply(.x, width = window_size, FUN = sd,   fill = NA, align = "right", partial = F),
      roll_acf1 = ~ rollapply(.x, width = window_size, FUN = acf1_fun, fill = NA, align = "right", partial = F),
      roll_cv = ~ rollapply(.x, widt = window_size, FUN = k_cv, fill = NA, align = "right", partial = F)
    ),
    .names = "{.col}_{.fn}"
  )) %>%
  ungroup() %>%
  select(year, Region, shannon_roll_mean:shannon_roll_cv)


rolling_metrics$Region <- factor(rolling_metrics$Region, levels=c('NFLS', 'SS', 'NEUS'))

gg_entropy_mean <- ggplot(rolling_metrics, aes(x = year, y = shannon_roll_mean, color = Region)) +
  facet_wrap(~ Region) +
  geom_line(linewidth = 2) +
  # geom_point(size = 3) +
  scale_color_viridis_d() +
  xlim(1970, 2005) +
  theme_bw(base_size = 14)

gg_entropy_mean

gg_entropy_sd <- ggplot(rolling_metrics, aes(x = year, y = shannon_roll_sd, color = Region)) +
  facet_wrap(~ Region) +
  geom_line(linewidth = 2) +
  # geom_point(size = 3) +
  scale_color_viridis_d() +
  xlim(1970, 2005) +
  theme_bw(base_size = 14)

gg_entropy_sd

gg_entropy_cv <- ggplot(rolling_metrics, aes(x = year, y = shannon_roll_cv, color = Region)) +
  geom_vline(data = df_gfb_windows, aes(xintercept = xintercept), linetype = "dashed", alpha = 0.8) +
  facet_wrap(~ Region, scale = "free") +
  geom_line(linewidth = 2) +
  # geom_point(size = 3) +
  labs(x = "Year",
       y = "Trophic \nStructural Stability (CV2)") +
  scale_color_viridis_d() +
  scale_x_continuous(limits = c(1970, 2005), breaks = scales::pretty_breaks()
  ) +
  scale_y_continuous(limits = c(0, 0.10), breaks = scales::pretty_breaks(), 
                     # expand = expansion(mult = c(0.25, 0.25))
  ) +
  theme_bw(base_size = 18) +
  theme(legend.position = "none")

gg_entropy_cv

# ggsave("../Scale_Dependence_Structural_Stability_LMEs/Figures/Figure SX - Windowed Shannon Structural and Stability.jpeg", plot = gg_entropy_cv, dpi = 300, width = 10, height = 4)

gg_entropy_acf1 <- ggplot(rolling_metrics, aes(x = year, y = shannon_roll_acf1, color = Region)) +
  facet_wrap(~ Region) +
  geom_line(linewidth = 2) +
  # geom_point(size = 3) +
  scale_color_viridis_d() +
  xlim(1970, 2005) +
  theme_bw(base_size = 14)

gg_entropy_acf1

gg_stability_grid <- plot_grid(gg_pred_precis, gg_entropy_cv,
                               ncol = 1, align = 'hv')
gg_stability_grid

# ggsave("../Scale_Dependence_Structural_Stability_LMEs/Figures/Figure SX - Windowed Structural Stability (Precision + CV2).jpeg", plot = gg_stability_grid, dpi = 300, width = 10, height = 9)


gg_entropy_grid <- plot_grid(gg_lme_entropy, gg_entropy_cv, 
                             ncol = 1, align = "hv")

gg_entropy_grid
# ggsave("../Scale_Dependence_Structural_Stability_LMEs/Figures/Figure SX - Windowed Shannon Structural Diversity and Stability.jpeg", plot = gg_entropy_grid, dpi = 300, width = 8, height = 7)
