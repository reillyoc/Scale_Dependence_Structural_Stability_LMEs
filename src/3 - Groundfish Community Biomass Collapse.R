# Visualize groundfish community biomass time series

# Author(s): Reilly O'Connor
# Version: 2025-06-25

# Load Pkgs
library(tidyverse)
library(vegan)
library(zoo)
library(reshape2)
library(cowplot)
library(ggrepel)
library(beepr)
library(easystats)

source("../Scale_Dependence_Structural_Stability_LMEs/src/0 - Functions.R")

# load data
df_nfls <- read.csv("../Scale_Dependence_Structural_Stability_LMEs/Data/Trawl Data/NFLS_NoPelagics_NW_Atlantic_trawlsurveydata_NFLS.csv", header = T)

df_ss <- read.csv("../Scale_Dependence_Structural_Stability_LMEs/Data/Trawl Data/SS_4VWXsurveydata_Fisher.csv", header = T)

df_neus <- read.csv("../Scale_Dependence_Structural_Stability_LMEs/Data/Trawl Data/NEUS_NoPelagics_NW_Atlantic_trawlsurveydataNumONLY.csv")

##### NFLS Groundfish Community Biomass #####
df_nfls_sub <-  df_nfls %>%
  #below filters to look at smaller cod management zones only
  # filter(lat <= 55.2) %>% #filter out sites North of 2J
  # filter(long <= -46.2) %>% #filter out sites East of 3L
  # mutate(keep = ifelse(lat < 48 & long < -54.5, FALSE, TRUE)) %>%
  # filter(keep == "TRUE") %>%
  # filter(lat >= 46.00) %>%g
  distinct(year, towid, spcode, weightpertow) %>%
  pivot_wider(names_from = spcode, values_from = weightpertow, values_fill = 0) %>%
  dplyr::select(year, towid, `73`, `74`, `75`) %>%
  pivot_longer(cols = c("73", "74", "75"), names_to = "spcode", values_to = "weightpertow") %>%
  arrange(year) %>%
  group_by(year, spcode) %>%
  summarize(geom_dens = CalcZeroInfGeomDens(weightpertow)) %>%
  group_by(spcode) %>%
  mutate(scale_geom_dens = scale(geom_dens),
         Region = "NFLS") %>%
  rename(Year = year) %>%
  filter(Year > 1972)

ggplot(data = df_nfls_sub, aes(x = Year, y = geom_dens, color = spcode)) +
  geom_vline(xintercept = 1989) +
  geom_vline(xintercept = 1993) +
  geom_point() +
  geom_line() +
  facet_wrap(~ spcode, scale = "free") +
  theme_classic() 

ggplot(data = df_nfls_sub, aes(x = Year, y = geom_dens, color = spcode)) +
  geom_vline(xintercept = 1989) +
  geom_vline(xintercept = 1993) +
  geom_point() +
  geom_line() +
  theme_classic()

df_nfls_sum <- df_nfls_sub %>%
  group_by(Year, Region) %>%
  summarise(sum_geom_dens = sum(geom_dens)) %>%
  ungroup() %>%
  mutate(scale_sum_geom_dens = scale(sum_geom_dens),
         spcode = "sum")

ggplot(data = df_nfls_sum, aes(x = Year, y = sum_geom_dens)) +
  geom_vline(xintercept = 1990) +
  geom_vline(xintercept = 1993) +
  geom_point() +
  geom_line() +
  theme_classic()

##### SS Groundfish Community Biomass #####
df_ss_sub <-  df_ss %>%
  distinct(year, towid, spcode, weightpertow) %>%
  pivot_wider(names_from = spcode, values_from = weightpertow, values_fill = 0) %>%
  dplyr::select(year, towid, `73`, `74`, `75`) %>%
  pivot_longer(cols = c("73", "74", "75"), names_to = "spcode", values_to = "weightpertow") %>%
  arrange(year) %>%
  group_by(year, spcode) %>%
  summarize(geom_dens = CalcZeroInfGeomDens(weightpertow)) %>%
  group_by(spcode) %>%
  mutate(scale_geom_dens = scale(geom_dens),
         Region = "SS") %>%
  rename(Year = year) %>%
  filter(Year > 1972 & Year < 2004)


ggplot(data = df_ss_sub, aes(x = Year, y = geom_dens, color = spcode)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ spcode, scale = "free") +
  theme_classic() 

ggplot(data = df_ss_sub, aes(x = Year, y = geom_dens, color = spcode)) +
  geom_point() +
  geom_line() +
  theme_classic()

df_ss_sum <- df_ss_sub %>%
  group_by(Year, Region) %>%
  summarise(sum_geom_dens = sum(geom_dens)) %>%
  ungroup() %>%
  mutate(scale_sum_geom_dens = scale(sum_geom_dens),
         spcode = "sum")

ggplot(data = df_ss_sum, aes(x = Year, y = sum_geom_dens)) +
  # geom_vline(xintercept = 1992) +
  geom_vline(xintercept = 1995) +
  geom_point() +
  geom_line() +
  theme_classic()


##### NEUS #####
df_neus_sub <- df_neus %>%
  filter(year >= 1973 & year < 2004 & nopertow > 0 & Trophic.level >= 3) %>%
  filter(! (MTW.Status == "Non-matching")) %>%
  mutate(
    nopertow = ifelse(nopertow == 0, 1, nopertow),
    towid = paste(lat, long, sep = "_")
  ) %>%
  filter(!is.na(nopertow)) %>%
  distinct(year, towid, spcode, .keep_all = TRUE) %>%
  pivot_wider(names_from = spcode, values_from = weightpertow, values_fill = 0) %>%
  dplyr::select(year, towid, `73`, `74`, `75`) %>%
  pivot_longer(cols = c("73", "74", "75"), names_to = "spcode", values_to = "weightpertow") %>%
  arrange(year) %>%
  group_by(year, spcode) %>%
  summarize(geom_dens = CalcZeroInfGeomDens(weightpertow)) %>%
  group_by(spcode) %>%
  mutate(scale_geom_dens = scale(geom_dens),
         Region = "NEUS") %>%
  rename(Year = year) 

ggplot(data = df_neus_sub, aes(x = Year, y = geom_dens, color = spcode)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ spcode, scale = "free") +
  theme_classic() 

ggplot(data = df_neus_sub, aes(x = Year, y = geom_dens, color = spcode)) +
  geom_point() +
  geom_line() +
  theme_classic()

df_neus_sum <- df_neus_sub %>%
  group_by(Year, Region) %>%
  summarise(sum_geom_dens = sum(geom_dens)) %>%
  ungroup() %>%
  mutate(scale_sum_geom_dens = scale(sum_geom_dens),
         spcode = "sum")

ggplot(data = df_neus_sum, aes(x = Year, y = log10(sum_geom_dens))) +
  geom_vline(xintercept = 1981) +
  geom_vline(xintercept = 1983) +
  geom_vline(xintercept = 1981) +
  geom_vline(xintercept = 1983) +
  geom_point() +
  geom_line() +
  theme_classic()

df_sum_gf <- rbind(df_nfls_sum, df_ss_sum, df_neus_sum)

#Regional Cod Spawning Stock Biomass Collapses
df_gfb_windows <- tribble(
  ~Region, ~xintercept,
  "NFLS",  1990,
  # "NFLS",  1994,
  "SS",    1992,
  # "SS",    1998,
  # "NEUS",  1981,
  "NEUS",  1990,
  # "NEUS",  1990,
  # "NEUS",  1993
)


df_sum_gf$Region <- factor(df_sum_gf$Region, levels=c('NFLS', 'SS', 'NEUS'))
df_gfb_windows$Region <- factor(df_gfb_windows$Region, levels=c('NFLS', 'SS', 'NEUS'))


gg_dens <- ggplot(data = df_sum_gf, aes(x = Year, y = log10(sum_geom_dens), color = Region)) +
  geom_vline(data = df_gfb_windows, aes(xintercept = xintercept), linetype = "dashed", alpha = 0.8) +
  # geom_vline(xintercept = 1981) +
  # geom_vline(xintercept = 1983) +
  facet_wrap(~ Region) +
  # geom_point() +
  scale_color_viridis_d() +
  geom_line(linewidth = 2) +
  scale_y_continuous(limits = c(-1.25, 1.75), breaks = scales::pretty_breaks()) +
  scale_x_continuous(limits = c(1970, 2005), breaks = scales::pretty_breaks()) +
  theme_bw(base_size = 14)

gg_dens

gg_dens_sc <- ggplot(data = df_sum_gf, aes(x = Year, y = log10(sum_geom_dens), color = Region)) +
  geom_vline(data = df_gfb_windows, aes(xintercept = xintercept), linetype = "dashed", alpha = 0.8) +
  # geom_vline(xintercept = 1981) +
  # geom_vline(xintercept = 1983) +
  facet_wrap(~ Region, scale = "free") +
  # geom_point() +
  scale_color_viridis_d() +
  geom_line(linewidth = 2) +
  # scale_y_continuous(limits = c(-1.25, 1.75), breaks = scales::pretty_breaks()) +
  scale_x_continuous(limits = c(1970, 2005), breaks = scales::pretty_breaks()) +
  theme_bw(base_size = 14)

gg_dens_sc

# ggsave("../Scale_Dependence_Structural_Stability_LMEs/Figures/Figure SX - Geometric Sum Mean Biomass Densities.jpeg", plot = gg_dens, width = 10, height = 3, dpi = 300)
# 
# ggsave("../Scale_Dependence_Structural_Stability_LMEs/Figures/Figure SX - Geometric Sum Mean Biomass Densities, free scale.jpeg", plot = gg_dens_sc, width = 10, height = 3, dpi = 300)

# write.csv(df_sum_gf, "../Scale_Dependence_Structural_Stability_LMEs/Data/Output Data/groundfish sum biomass.csv")

