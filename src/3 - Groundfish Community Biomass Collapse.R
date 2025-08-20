# NMDS on Annual Species Assemblages from 1981 onwards...

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

source("../Structural_Stability_LMEs/src/Functions.R")

# load data
df_nfls <- read.csv("../Structural_Stability_LMEs/data/NFLS_NoPelagics_NW_Atlantic_trawlsurveydata.csv", header = T)

df_ss <- read.csv("../Structural_Stability_LMEs/data/SS_4VWXsurveydata_Fisher.csv", header = T)

df_neus <- read.csv("../Structural_Stability_LMEs/data/NEUS_NoPelagics_NW_Atlantic_trawlsurveydataNumONLY.csv")


##### NFLS Groundfish Community Biomass #####
df_nfls_sub <-  df_nfls %>%
  filter(lat <= 55.2) %>% # filter out sites North of 2J
  filter(long <= -46.2) %>% # filter out sites East of 3L
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

library(nlme)

lm_nfls_sum <- gls(sum_geom_dens ~ Year, data = df_nfls_sum,
                   correlation = corAR1(form =~ Year))
summary(lm_nfls_sum)

plot(df_nfls_sum$resid_sum_dens)

acf(df_nfls_sum$sum_geom_dens)

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
  geom_vline(xintercept = 1992) +
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
  "NFLS",  1993,
  "SS",    1992,
  "SS",    1998,
  "NEUS",  1981,
  "NEUS",  1991,
  # "NEUS",  1990,
  # "NEUS",  1993
)


df_sum_gf$Region <- factor(df_sum_gf$Region, levels=c('NFLS', 'SS', 'NEUS'))
df_gfb_windows$Region <- factor(df_gfb_windows$Region, levels=c('NFLS', 'SS', 'NEUS'))


gg_dens <- ggplot(data = df_sum_gf, aes(x = Year, y = log10(sum_geom_dens), color = Region)) +
  geom_vline(data = df_gfb_windows, aes(xintercept = xintercept), linetype = "dashed", alpha = 0.8) +
  # geom_vline(xintercept = 1981) +
  # geom_vline(xintercept = 1983) +
  facet_wrap(~ Region, scale = "free") +
  # geom_point() +
  scale_color_viridis_d() +
  geom_line(linewidth = 2) +
  xlim(1970, 2006) +
  theme_bw(base_size = 14)

gg_dens

ggsave("../Structural_Stability_LMEs/figures/Figure SX - Geometric Sum Mean Biomass Densities.jpeg", plot = gg_dens, width = 10, height = 3, dpi = 300)


cowplot::plot_grid(gg_dens, gg_pred_precis, nrow = 2)

# write.csv(df_sum_gf, "../Structural_Stability_LMEs/data/groundfish sum biomass.csv")

