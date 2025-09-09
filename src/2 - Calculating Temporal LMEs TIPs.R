# Structural Stability in Large Marine Ecosystems

# Author(s): Reilly O'Connor
# Version: 2025-05-07

# Load Pkgs
library(tidyverse)
library(vegan)
library(zoo)
library(reshape2)
library(cowplot)
library(beepr)
library(easystats)

source("../Scale_Dependence_Structural_Stability_LMEs/src/Functions.R")

# load data
df_nfls <- read.csv("../Scale_Dependence_Structural_Stability_LMEs/Data/Trawl Data/NFLS_NoPelagics_NW_Atlantic_trawlsurveydata_NFLS.csv", header = T)

df_ss <- read.csv("../Scale_Dependence_Structural_Stability_LMEs/Data/Trawl Data/SS_4VWXsurveydata_Fisher.csv", header = T)

df_neus <- read.csv("../Scale_Dependence_Structural_Stability_LMEs/Data/Trawl Data/NEUS_NoPelagics_NW_Atlantic_trawlsurveydataNumONLY.csv")


##### Species Accumulation Curves by Year #####
#NFLS
df_sp_accum_by_year_nfls <- df_nfls %>%
  filter(lat <= 55.2) %>% # filter out sites North of 2J
  filter(long <= -46.2) %>% # filter out sites East of 3L
  # mutate(keep = ifelse(lat < 48 & long < -54.5, FALSE, TRUE)) %>%
  # filter(keep == "TRUE") %>%
  # filter(lat >= 46.00) %>%
  dplyr::select(towid, year, spcode, nopertow) %>%
  filter(year >= 1973) %>%
  group_by(towid, year, spcode) %>%
  reframe(pres = as.integer(mean(nopertow) > 0)) %>%  # <-- Convert to presence/absence
  ungroup() %>%
  group_by(year) %>%
  nest() %>%
  mutate(
    # pivot each year’s data to wide (sites × species), filling NAs with 0
    sp_mat = map(data, ~ 
                   pivot_wider(.x, names_from = spcode, values_from = pres, values_fill = 0) %>% 
                   dplyr::select(-towid)
    ),
    # run the random accumulation curve on presence/absence matrix
    sac = map(sp_mat, ~ specaccum(.x, method = "random")),
    
    # pull out the sites, mean richness and sd into a tibble
    sac_df = map(sac, ~ tibble(
      sites    = .$sites,
      richness = .$richness,
      sd       = .$sd
    ))
  ) %>%
  dplyr::select(year, sac_df) %>%
  unnest(sac_df)


gg_nfls_sac <- ggplot(df_sp_accum_by_year_nfls, aes(x = sites, y = richness)) +
  geom_ribbon(aes(ymin = richness - sd, ymax = richness + sd), fill = "steelblue") +
  geom_line(size = 1) +
  facet_wrap(~ year, nrow = 4) +
  labs(x = "Number of Tows (sites)",
       y = "Species Richness",
       title = "Species–accumulation curves by year") +
  theme_classic()

gg_nfls_sac

df_sites_by_year_nfls <- df_sp_accum_by_year_nfls %>%
  dplyr::select(year, sites) %>%
  group_by(year) %>%
  reframe(max_sites = max(sites))

ggplot(df_sites_by_year_nfls, aes(x = year, y = max_sites)) +
  geom_point() +
  xlim(1981, 2004) +
  geom_line(size = 1) +
  labs(x = "Year",
       y = "Number of Sites") +
  theme_classic() 

#Scotian Shelf

#Transpose Data for use with Vegan specaccum()
df_sp_accum_by_year <- df_ss %>%
  dplyr::select(towid, year, spcode, nopertow) %>%
  filter(year >= 1973 & year < 2004) %>%
  group_by(year) %>%
  nest() %>% 
  mutate(
    # pivot each year’s data to wide (sites × species), filling NAs with 0
    sp_mat = map(data, ~ 
                   pivot_wider(.x, names_from = spcode, values_from = nopertow, values_fill = 0) %>% 
                   dplyr::select(-towid)
    ),
    # run the random accumulation curve
    sac   = map(sp_mat, ~ specaccum(.x, method = "random")),
    # pull out the sites, mean richness and sd into a tibble
    sac_df = map(sac, ~ tibble(
      sites    = .$sites,
      richness = .$richness,
      sd       = .$sd
    ))
  ) %>%
  dplyr::select(year, sac_df) %>%
  unnest(sac_df)

gg_ss_sac <- ggplot(df_sp_accum_by_year, aes(x = sites, y = richness)) +
  geom_ribbon(aes(ymin = richness - sd, ymax = richness + sd), fill = "steelblue") +
  geom_line(size = 1) +
  facet_wrap(~ year, nrow = 4) +
  labs(x = "Number of Tows (sites)",
       y = "Species Richness",
       title = "Species–accumulation curves by year") +
  theme_classic()

gg_ss_sac

df_sites_by_year <- df_sp_accum_by_year %>%
  dplyr::select(year, sites) %>%
  group_by(year) %>%
  reframe(max_sites = max(sites))

ggplot(df_sites_by_year, aes(x = year, y = max_sites)) +
  geom_point() +
  geom_line(size = 1) +
  labs(x = "Year",
       y = "Number of Sites") +
  theme_classic() 

#NEUS
#Transpose Data for use with Vegan specaccum()
df_sp_accum_by_year_neus <- df_neus %>%
  dplyr::select(year, lat, long, spcode, nopertow) %>%
  filter(year >= 1973 & year < 2004) %>%
  mutate(
    nopertow = ifelse(nopertow == 0, 1, nopertow),
    towid = paste(lat, long, sep = "_")
  ) %>%
  filter(!is.na(nopertow)) %>%
  distinct(year, towid, spcode, .keep_all = TRUE) %>%  # Ensure uniqueness
  dplyr::select(year, towid, spcode, nopertow) %>%
  group_by(year) %>%
  nest() %>% 
  mutate(
    sp_mat = map(data, ~ 
                   pivot_wider(
                     .x,
                     names_from = spcode,
                     values_from = nopertow,
                     values_fill = list(nopertow = 0)
                   ) %>% 
                   dplyr::select(-towid)
    ),
    sac = map(sp_mat, ~ specaccum(.x, method = "random")),
    sac_df = map(sac, ~ tibble(
      sites    = .$sites,
      richness = .$richness,
      sd       = .$sd
    ))
  ) %>%
  dplyr::select(year, sac_df) %>%
  unnest(sac_df)

gg_neus_sac <- ggplot(df_sp_accum_by_year_neus, aes(x = sites, y = richness)) +
  geom_ribbon(aes(ymin = richness - sd, ymax = richness + sd), fill = "steelblue") +
  geom_line(size = 1) +
  facet_wrap(~ year, nrow = 4) +
  labs(x = "Number of Tows (sites)",
       y = "Species Richness",
       title = "Species–accumulation curves by year") +
  theme_classic()

gg_neus_sac

df_sites_by_year_neus <- df_sp_accum_by_year_neus %>%
  dplyr::select(year, sites) %>%
  group_by(year) %>%
  reframe(max_sites = max(sites))

ggplot(df_sites_by_year_neus, aes(x = year, y = max_sites)) +
  geom_point() +
  geom_line(size = 1) +
  labs(x = "Year",
       y = "Number of Sites") +
  theme_classic() 

##### Total Species Richness per Year #####
df_nfls_sp <- df_sp_accum_by_year_nfls %>%
  group_by(year) %>%
  filter(sites == max(sites)) %>%
  rename(total_sr = richness) %>%
  dplyr::select(-sd)

df_ss_sp <- df_sp_accum_by_year %>%
  group_by(year) %>%
  filter(sites == max(sites)) %>%
  rename(total_sr = richness) %>%
  dplyr::select(-sd)

df_neus_sp <- df_sp_accum_by_year_neus %>%
  group_by(year) %>%
  filter(sites == max(sites)) %>%
  rename(total_sr = richness) %>%
  dplyr::select(-sd)

ggplot() +
  geom_vline(xintercept = 1990, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 1995, linetype = "dashed", alpha = 0.5) +
  xlim(1977, 2003) +
  geom_point(data = df_nfls_sp, aes(x = year, y = scale(total_sr))) +
  geom_line(data = df_nfls_sp, aes(x = year, y = scale(total_sr)), color = "navy") +
  geom_point(data = df_nfls_sp, aes(x = year, y = scale(sites))) +
  geom_line(data = df_nfls_sp, aes(x = year, y = scale(sites)), color = "orange2") +
  labs(x = "Year",
       y = "Total Species Richness") +
  #xlim(1981, 2004) +
  theme_classic()


ggplot() +
  geom_vline(xintercept = 1995, linetype = "dashed", alpha = 0.5) +
  geom_point(data = df_ss_sp, aes(x = year, y = scale(total_sr))) +
  geom_line(data = df_ss_sp, aes(x = year, y = scale(total_sr)), color = "navy") +
  geom_point(data = df_ss_sp, aes(x = year, y = scale(sites))) +
  geom_line(data = df_ss_sp, aes(x = year, y = scale(sites)), color = "orange2") +
  labs(x = "Year",
       y = "Total Species Richness") +
  theme_classic()


ggplot() +
  geom_vline(xintercept = 1995, linetype = "dashed", alpha = 0.5) +
  geom_point(data = df_neus_sp, aes(x = year, y = scale(total_sr))) +
  geom_line(data = df_neus_sp, aes(x = year, y = scale(total_sr)), color = "navy") +
  geom_point(data = df_neus_sp, aes(x = year, y = scale(sites))) +
  geom_line(data = df_neus_sp, aes(x = year, y = scale(sites)), color = "orange2") +
  labs(x = "Year",
       y = "Total Species Richness") +
  theme_classic()

df_nfls_sp$Region <- "NFLS"
df_ss_sp$Region <- "SS"
df_neus_sp$Region <- "NEUS"

df_lme_sp <- rbind(df_nfls_sp, df_ss_sp, df_neus_sp)

# write.csv(df_lme_sp, "../Scale_Dependence_Structural_Stability_LMEs/Data/Output Data/number_sites_regions.csv")

ggsave("../Scale_Dependence_Structural_Stability_LMEs/Figures/Figure S2 - NFLS Species Accumulation Curves.jpeg", plot = gg_nfls_sac, dpi = 300, width = 14, height = 7)

ggsave("../Scale_Dependence_Structural_Stability_LMEs/Figures/Figure S3 - SS Species Accumulation Curves.jpeg", plot = gg_ss_sac, dpi = 300, width = 14, height = 7)

ggsave("../Scale_Dependence_Structural_Stability_LMEs/Figures/Figure S4 - NEUS Species Accumulation Curves.jpeg", plot = gg_neus_sac, dpi = 300, width = 14, height = 7)

##### Quantifying Species Richness by Trophic Interval & TIPs #####
#NFLS
df_nfls_pa <- df_nfls  %>%
  # filter(lat <= 55.2) %>% # filter out sites North of 2J 
  # filter(long <= -46.2) %>% # filter out sites East of 3L
  # mutate(keep = ifelse(lat < 48 & long < -54.5, FALSE, TRUE)) %>%
  # filter(keep == "TRUE") %>%
  #filter(lat >= 46.00) %>%
  filter(year >= 1973 & nopertow > 0) %>%
  mutate(TL = round(Trophic.level, 1),
         Trophic_interval = case_when(
           TL >= 3.0  & TL < 3.5  ~ "3.0-3.4",
           TL >= 3.5  & TL < 4.0  ~ "3.5-3.9",
           TL >= 4.0  & TL < 4.5  ~ "4.0-4.4",
           TL >= 4.5  & TL <= 4.9  ~ "4.5-4.9",
           TRUE ~ NA_character_)) %>%
  group_by(towid, year, Trophic_interval, spcode) %>%
  summarize(pres = as.integer(mean(nopertow) > 0), .groups = "drop")

df_nfls_chao2_ti <- tibble()

Years <- df_nfls_pa %>%
  distinct(year) %>%
  arrange(year) %>%
  pull()

trophic <- df_nfls_pa %>%
  distinct(Trophic_interval) %>%
  pull()

#Group Species by Trophic Interval - Estimate Spp Richness by Trophic Interval
for (ti in trophic) {
  for (time in Years) {
    
    df_nfls_pa_sub <- df_nfls_pa %>%
      filter(year == time & Trophic_interval == ti) %>%
      dplyr::select(-year, -Trophic_interval) %>%
      pivot_wider(names_from = spcode, values_from = pres, values_fill = 0) %>%
      dplyr::select(-towid)
    
    if (nrow(df_nfls_pa_sub) < 2 || ncol(df_nfls_pa_sub) == 0) next  # skip empty
    
    # --- Chao2 estimates ---
    temp_chao2 <- specpool(df_nfls_pa_sub) %>%
      as_tibble() %>%
      mutate(year = time,
             Trophic_interval = ti)
    
    df_nfls_chao2_ti <- bind_rows(df_nfls_chao2_ti, temp_chao2)
    
  }
}


df_nfls_sp_tips <- df_nfls_chao2_ti %>%
  complete(
    year = full_seq(df_nfls$year[df_nfls$year >= 1973], 1),
    Trophic_interval = c("3.0-3.4", "3.5-3.9", "4.0-4.4", "4.5-4.9"),
    fill = list(Species = 0, chao = 0, jack1 = 0, jack2 = 0, boot = 0)
  )  %>%
  group_by(year) %>%
  mutate(
    total_sr = sum(Species),
    total_sr_chao = sum(chao),
    total_sr_jack1 = sum(jack1),
    total_sr_jack2 = sum(jack2),
    total_sr_boot = sum(boot),
    tips = Species/total_sr,
    tips_chao = chao/total_sr_chao,
    tips_jack1 = jack1/total_sr_jack1,
    tips_jack2 = jack2/total_sr_jack2,
    tips_boot = boot/total_sr_boot
  ) %>%
  ungroup()


ggplot(df_nfls_sp_tips, aes(x = year, y = total_sr)) +
  geom_point() +
  geom_line() +
  labs(x = "Year",
       y = "Total Species Richness") +
  theme_classic()

ggplot(df_nfls_sp_tips, aes(x = year, y = Species, color = Trophic_interval)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("#4245c4", "#f78c2f",  "#649f4d", "#b61790")) +
  labs(x = "Year",
       y = "Species Richness") +
  #xlim(1981, 2004) +
  theme_classic()

ggplot(df_nfls_sp_tips, aes(x = year, y = tips_chao, color = Trophic_interval)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("#4245c4", "#f78c2f",  "#649f4d", "#b61790")) +
  labs(x = "Year",
       y = "Trophic Interval Proportions") +
  #xlim(1981, 2004) +
  theme_classic()

#SS
df_ss_pa <- df_ss %>%
  filter(year >= 1973 & year < 2004 & nopertow > 0) %>%
  mutate(TL = round(Trophic.level, 1),
         Trophic_interval = case_when(
           TL >= 3.0  & TL < 3.5  ~ "3.0-3.4",
           TL >= 3.5  & TL < 4.0  ~ "3.5-3.9",
           TL >= 4.0  & TL < 4.5  ~ "4.0-4.4",
           TL >= 4.5  & TL <= 4.9  ~ "4.5-4.9",
           TRUE ~ NA_character_)) %>%
  group_by(towid, year, Trophic_interval, spcode) %>%
  summarize(pres = as.integer(mean(nopertow) > 0), .groups = "drop")

df_ss_chao2_ti <- tibble()

Years <- df_ss_pa %>%
  distinct(year) %>%
  arrange(year) %>%
  pull()

trophic <- df_ss_pa %>%
  distinct(Trophic_interval) %>%
  pull()

#Group Species by Trophic Interval - Estimate Spp Richness by Trophic Interval
for (ti in trophic) {
  for (time in Years) {
    
    df_ss_pa_sub <- df_ss_pa %>%
      filter(year == time & Trophic_interval == ti) %>%
      dplyr::select(-year, -Trophic_interval) %>%
      pivot_wider(names_from = spcode, values_from = pres, values_fill = 0) %>%
      dplyr::select(-towid)
    
    if (nrow(df_ss_pa_sub) < 2 || ncol(df_ss_pa_sub) == 0) next  # skip empty
    
    # --- Chao2 estimates ---
    temp_chao2 <- specpool(df_ss_pa_sub) %>%
      as_tibble() %>%
      mutate(year = time,
             Trophic_interval = ti)
    
    df_ss_chao2_ti <- bind_rows(df_ss_chao2_ti, temp_chao2)
    
  }
}


df_ss_sp_tips <- df_ss_chao2_ti %>%
  group_by(year) %>%
  mutate(
    total_sr = sum(Species),
    total_sr_chao = sum(chao),
    total_sr_jack1 = sum(jack1),
    total_sr_jack2 = sum(jack2),
    total_sr_boot = sum(boot),
    tips = Species/total_sr,
    tips_chao = chao/total_sr_chao,
    tips_jack1 = jack1/total_sr_jack1,
    tips_jack2 = jack2/total_sr_jack2,
    tips_boot = boot/total_sr_boot
  ) %>%
  ungroup()

ggplot(df_ss_sp_tips, aes(x = year, y = total_sr)) +
  geom_point() +
  geom_line() +
  labs(x = "Year",
       y = "Total Species Richness") +
  theme_classic()

ggplot(df_ss_sp_tips, aes(x = year, y = Species, color = Trophic_interval)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("#4245c4", "#f78c2f",  "#649f4d", "#b61790")) +
  labs(x = "Year",
       y = "Species Richness") +
  theme_classic()

ggplot(df_ss_sp_tips, aes(x = year, y = tips_chao, color = Trophic_interval)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("#4245c4", "#f78c2f",  "#649f4d", "#b61790")) +
  labs(x = "Year",
       y = "Trophic Interval Proportions") +
  theme_classic()

#NEUS
df_neus_pa <- df_neus  %>%
  filter(year >= 1973 & year < 2004 & nopertow > 0 & Trophic.level >= 3) %>%
  filter(! (MTW.Status == "Non-matching")) %>%
  mutate(
    nopertow = ifelse(nopertow == 0, 1, nopertow),
    towid = paste(lat, long, sep = "_")
  ) %>%
  filter(!is.na(nopertow)) %>%
  distinct(year, towid, spcode, .keep_all = TRUE) %>%  
  dplyr::select(year, towid, spcode, nopertow, Trophic.level) %>%
  mutate(TL = round(Trophic.level, 1),
         Trophic_interval = case_when(
           TL >= 3.0  & TL < 3.5  ~ "3.0-3.4",
           TL >= 3.5  & TL < 4.0  ~ "3.5-3.9",
           TL >= 4.0  & TL < 4.5  ~ "4.0-4.4",
           TL >= 4.5  & TL <= 4.9  ~ "4.5-4.9",
           TRUE ~ NA_character_)) %>%
  group_by(towid, year, Trophic_interval, spcode) %>%
  summarize(pres = as.integer(mean(nopertow) > 0), .groups = "drop")

df_neus_chao2_ti <- tibble()

Years <- df_neus_pa %>%
  distinct(year) %>%
  arrange(year) %>%
  pull()

trophic <- df_neus_pa %>%
  distinct(Trophic_interval) %>%
  pull()

#Group Species by Trophic Interval - Estimate Spp Richness by Trophic Interval
for (ti in trophic) {
  for (time in Years) {
    
    df_neus_pa_sub <- df_neus_pa %>%
      filter(year == time & Trophic_interval == ti) %>%
      dplyr::select(-year, -Trophic_interval) %>%
      pivot_wider(names_from = spcode, values_from = pres, values_fill = 0) %>%
      dplyr::select(-towid)
    
    if (nrow(df_neus_pa_sub) < 2 || ncol(df_neus_pa_sub) == 0) next  # skip empty
    
    # --- Chao2 estimates ---
    temp_chao2 <- specpool(df_neus_pa_sub) %>%
      as_tibble() %>%
      mutate(year = time,
             Trophic_interval = ti)
    
    df_neus_chao2_ti <- bind_rows(df_neus_chao2_ti, temp_chao2)
    
  }
}


df_neus_sp_tips <- df_neus_chao2_ti %>%
  group_by(year) %>%
  mutate(
    total_sr = sum(Species),
    total_sr_chao = sum(chao),
    total_sr_jack1 = sum(jack1),
    total_sr_jack2 = sum(jack2),
    total_sr_boot = sum(boot),
    tips = Species/total_sr,
    tips_chao = chao/total_sr_chao,
    tips_jack1 = jack1/total_sr_jack1,
    tips_jack2 = jack2/total_sr_jack2,
    tips_boot = boot/total_sr_boot
  ) %>%
  ungroup()


ggplot(df_neus_sp_tips, aes(x = year, y = total_sr)) +
  geom_point() +
  geom_line() +  
  labs(x = "Year",
       y = "Total Species Richness") +
  
  theme_classic()

ggplot(df_neus_sp_tips, aes(x = year, y = Species, color = Trophic_interval)) +
  geom_point() +
  geom_line() +  
  scale_color_manual(values = c("#4245c4", "#f78c2f",  "#649f4d", "#b61790")) +
  labs(x = "Year",
       y = "Species Richness") +
  
  theme_classic()

ggplot(df_neus_sp_tips, aes(x = year, y = tips, color = Trophic_interval)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("#4245c4", "#f78c2f",  "#649f4d", "#b61790")) +
  labs(x = "Year",
       y = "Trophic Interval Proportions") +
  theme_classic()



df_nfls_sp_tips$Region <- "NFLS"
df_ss_sp_tips$Region <- "SS"
df_neus_sp_tips$Region <- "NEUS"

df_sp_tips <- rbind(df_nfls_sp_tips, df_ss_sp_tips, df_neus_sp_tips)


#ACF Plots for TIPs
df_nfls_sp_tips_sub <- df_sp_tips %>%
  filter(Region == "NFLS") %>%
  dplyr::select(year, Region, Trophic_interval, tips_chao) %>%
  rename(tips = tips_chao)

df_sp_tips_nss <- df_sp_tips %>%
  filter(Region == "SS" | Region == "NEUS") %>%
  dplyr::select(year, Region, Trophic_interval, tips)

df_sp_tips_sub <- rbind(df_nfls_sp_tips_sub, df_sp_tips_nss)

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


acf_df$Region <- factor(acf_df$Region, levels=c('NFLS', 'SS', 'NEUS'))
ci_df$Region <- factor(ci_df$Region, levels=c('NFLS', 'SS', 'NEUS'))

gg_tips_acf <- ggplot(acf_df, aes(x = as.factor(lag), y = acf)) +  
  geom_hline(data = ci_df, aes(yintercept = ci), linetype = "dashed", color = "blue", alpha = 0.5) +
  geom_hline(data = ci_df, aes(yintercept = -ci), linetype = "dashed", color = "blue", alpha = 0.5) +
  geom_hline(yintercept = 0, color = "black", alpha = 0.5) +
  geom_col() +
  facet_wrap(~Region + Trophic_interval) +
  theme_bw(base_size = 14) +
  labs(x = "Lag", y = "ACF")

gg_tips_acf

ggsave("../Scale_Dependence_Structural_Stability_LMEs/Figures/Figure SX - TIPs ACF.jpeg", plot = gg_tips_acf, dpi = 300, width = 10, height = 9)

write.csv(df_sp_tips, "../Scale_Dependence_Structural_Stability_LMEs/Data/Output Data/richness_tips_split_chao.csv")

#Look at observed TIPs

df_sp_tips$Region <- factor(df_sp_tips$Region, levels=c('NFLS', 'SS', 'NEUS'))

gg_tips <- ggplot(df_sp_tips, aes(x = year, y = tips, color = Trophic_interval)) +
  # geom_vline(xintercept = 1990, linetype = "dashed", alpha = 0.5) +
  # geom_vline(data = df_sp_tips %>% filter(Region == "NFLS" & year == 1995),
  # aes(xintercept = year),
  # linetype = "dashed",
  # alpha = 0.5) +
  #geom_point(size = 3, shape = 21, color = "black",  stroke = 0.5, alpha = 0.8, aes(fill = Trophic_interval)) +
  geom_line(alpha = 0.8, linewidth = 2) + 
  scale_color_manual(values = c("#4245c4", "#f78c2f",  "#649f4d", "#b61790")) +
  scale_fill_manual(values = c("#4245c4", "#f78c2f",  "#649f4d", "#b61790")) +
  labs(x = "Year",
       y = "Trophic Interval Proportions (TIPs)") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none") +
  ylim(-0.05, 0.65) +
  xlim(1973, 2004) +
  facet_wrap(~ Region) +
  theme(text = element_text(family = "Arial"), 
        legend.position = "none") 

gg_tips

# ggsave("../Scale_Dependence_Structural_Stability_LMEs/Figures/Figure 3 - Oberved Temporal TIPs.jpeg", plot = gg_tips, dpi = 300, width = 10, height = 4)


