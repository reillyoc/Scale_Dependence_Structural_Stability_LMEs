# Bray-Curtis Similarity and Community NMDS on Annual Species Assemblages

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
df_nfls_nm <- read.csv("../Scale_Dependence_Structural_Stability_LMEs/Data/Trawl Data/NFLS_NoPelagics_NW_Atlantic_trawlsurveydata_NFLS.csv", header = T)

df_ss_nm <- read.csv("../Scale_Dependence_Structural_Stability_LMEs/Data/Trawl Data/SS_4VWXsurveydata_Fisher.csv", header = T)

df_neus_nm <- read.csv("../Scale_Dependence_Structural_Stability_LMEs/Data/Trawl Data/NEUS_NoPelagics_NW_Atlantic_trawlsurveydataNumONLY.csv")

##### NFLS #####
df_nfls_sub <-  df_nfls_nm %>%
  filter(lat <= 55.2) %>% # filter out sites North of 2J
  filter(long <= -46.2) %>% # filter out sites East of 3L
  # mutate(keep = ifelse(lat < 48 & long < -54.5, FALSE, TRUE)) %>%
  # filter(keep == "TRUE") %>%
  # filter(lat >= 46.00) %>%
  dplyr::select(towid, year, spcode, nopertow) %>%
  filter(year >= 1973) %>%
  group_by(towid, year, spcode) %>%
  reframe(pres = as.integer(mean(nopertow) > 0)) %>%
  distinct(year, spcode, pres) %>%
  arrange(year) %>%
  pivot_wider(names_from = spcode, values_from = pres, values_fill = 0)

nfls_species_matrix <- df_nfls_sub %>% dplyr::select(-year)
years <- df_nfls_sub$year

rownames(nfls_species_matrix) <- years  

nfls_bray <- vegdist(nfls_species_matrix, method = "bray")

nfls_bray_mat <- as.matrix(nfls_bray)

nfls_bray_1981 <- tibble(Year = 1973:2004)
nfls_bray_1981$BC_dis <- nfls_bray_mat[,9]
nfls_bray_1981 <- nfls_bray_1981 %>%
  mutate(BC_simm = 1 - BC_dis,
         Region = "NFLS")

gg_nfls_bc <- ggplot(data = nfls_bray_1981, aes(x = Year, y = BC_simm)) +
  geom_point(size = 4, shape = 21, stroke = 0.5) +
  geom_line(linewidth = 2) +
  labs(title = "NFLS") +
  theme_bw(base_size = 14) +
  ylim(0.65, 1) +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        text = element_text(family = "Arial"), 
        legend.position = "none") 

gg_nfls_bc

nfls_nmds_result <- metaMDS(nfls_species_matrix, distance = "bray", k = 2, trymax = 100)

nfls_nmds_scores <- scores(nfls_nmds_result, display = "sites") %>% 
  as.data.frame() %>%
  mutate(Year = df_nfls_sub$year,
         Region = "NFLS") 

gg_nfls_nmds <- ggplot(data = nfls_nmds_scores, aes(x = NMDS1, y = NMDS2, label = Year)) +
  geom_path(color = "gray") +  # show temporal trajectory
  geom_point(size = 4, shape = 21, stroke = 0.5, fill = "aliceblue") +
  scale_fill_manual(values = c("#A6CAEC", "#62778B")) +  
  geom_text_repel(size = 4, max.overlaps = Inf) +  
  labs(title = "NFLS") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        text = element_text(family = "Arial"), 
        legend.position = "none") 

gg_nfls_nmds


##### SS #####
df_ss_sub <- df_ss_nm %>%
  filter(year >= 1973 & year < 2004 & nopertow > 0) %>%
  group_by(towid, year, spcode) %>%
  summarize(pres = as.integer(mean(nopertow) > 0), .groups = "drop") %>%
  distinct(year, spcode, pres) %>%
  arrange(year) %>%
  pivot_wider(names_from = spcode, values_from = pres, values_fill = 0)


ss_species_matrix <- df_ss_sub %>% dplyr::select(-year)

years <- df_ss_sub$year

rownames(ss_species_matrix) <- years

ss_bray <- vegdist(ss_species_matrix, method = "bray")

ss_bray_mat <- as.matrix(ss_bray)

ss_bray_1981 <- tibble(Year = 1973:2003)
ss_bray_1981$BC_dis <- ss_bray_mat[,9]
ss_bray_1981 <- ss_bray_1981 %>%
  mutate(BC_simm = 1 - BC_dis,
         Region = "SS")

gg_ss_bc <- ggplot(data = ss_bray_1981, aes(x = Year, y = BC_simm)) +
  geom_point(size = 4, shape = 21, stroke = 0.5) +
  geom_line() +
  labs(title = "SS") +
  theme_bw(base_size = 14) +
  ylim(0.65, 1) +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        text = element_text(family = "Arial"), 
        legend.position = "none") 

gg_ss_bc


ss_nmds_result <- metaMDS(ss_species_matrix, distance = "bray", k = 2, trymax = 100)

ss_nmds_scores <- scores(ss_nmds_result, display = "sites") %>% 
  as.data.frame() %>%
  mutate(Year = df_ss_sub$year,
         Region = "SS")

gg_ss_nmds <- ggplot(data = ss_nmds_scores, aes(x = NMDS1, y = NMDS2, label = Year)) +
  geom_path(color = "gray") +  # show temporal trajectory
  geom_point(size = 4, shape = 21, stroke = 0.5, fill = "aliceblue") +
  scale_fill_manual(values = c("#A6CAEC", "#62778B")) +  
  geom_text_repel(size = 4, max.overlaps = Inf) +  
  labs(title = "SS") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        text = element_text(family = "Arial"), 
        legend.position = "none") 

gg_ss_nmds



##### NEUS #####
df_neus_sub <- df_neus_nm %>%
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
  group_by(towid, year, spcode) %>%
  summarize(pres = as.integer(mean(nopertow) > 0), .groups = "drop") %>%
  distinct(year, spcode, pres) %>%
  arrange(year) %>%
  pivot_wider(names_from = spcode, values_from = pres, values_fill = 0)

neus_species_matrix <- df_neus_sub %>% dplyr::select(-year)

years <- df_neus_sub$year

rownames(neus_species_matrix) <- years

neus_bray <- vegdist(neus_species_matrix, method = "bray")

neus_bray_mat <- as.matrix(neus_bray)

neus_bray_1981 <- tibble(Year = 1973:2003)
neus_bray_1981$BC_dis <- neus_bray_mat[,9]
neus_bray_1981 <- neus_bray_1981 %>%
  mutate(BC_simm = 1 - BC_dis,
         Region = "NEUS")

gg_neus_bc <- ggplot(data = neus_bray_1981, aes(x = Year, y = BC_simm)) +
  geom_point(size = 4, shape = 21, stroke = 0.5) +
  geom_line() +
  labs(title = "SS") +
  theme_bw(base_size = 14) +
  ylim(0.6, 1) +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        text = element_text(family = "Arial"), 
        legend.position = "none") 

gg_neus_bc

neus_nmds_result <- metaMDS(neus_species_matrix, distance = "bray", k = 2, trymax = 100)

neus_nmds_scores <- scores(neus_nmds_result, display = "sites") %>% 
  as.data.frame() %>%
  mutate(Year = df_neus_sub$year,
         Region = "NEUS")

gg_neus_nmds <- ggplot(data = neus_nmds_scores, aes(x = NMDS1, y = NMDS2, label = Year)) +
  geom_path(color = "gray", aes(group = NA)) +  # show temporal trajectory
  geom_point(size = 4, shape = 21, stroke = 0.5, fill = "aliceblue") +
  geom_text_repel(size = 4, max.overlaps = Inf) +  
  labs(title = "NEUS") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        text = element_text(family = "Arial"), 
        legend.position = "none") 

gg_neus_nmds

gg_nmds <- plot_grid(gg_nfls_nmds, gg_ss_nmds, gg_neus_nmds, nrow = 1, align = "hv")

df_nmds_scores <- rbind(nfls_nmds_scores, ss_nmds_scores, neus_nmds_scores)

df_nmds_scores$Region <- factor(df_nmds_scores$Region, levels=c('NFLS', 'SS', 'NEUS'))

df_nmds_scores$Year <- as.integer(df_nmds_scores$Year)
df_nmds_scores_labels <- df_nmds_scores %>% filter(Year %% 3 == 0)


gg_nmds_all <-  ggplot(data = df_nmds_scores, aes(x = NMDS1, y = NMDS2, label = Year)) +
  geom_path(color = "gray", aes(group = NA)) +  # show temporal trajectory
  geom_point(size = 4, shape = 21, stroke = 0.5, fill = "aliceblue") +
  geom_text_repel(
    data = df_nmds_scores_labels,
  point.padding = 0.3,      
  box.padding = 0.5,        
  max.overlaps = Inf,
  segment.size = 0.3,     
  segment.color = "gray50", # more subtle connecting lines
  min.segment.length = 0.1) +
  labs(x = "NMDS 1",
       y = "NMDS 2") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none") +
  facet_wrap(~ Region) +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        text = element_text(family = "Arial"), 
        legend.position = "none") 

gg_nmds_all



df_bc_all <- rbind(nfls_bray_1981, ss_bray_1981, neus_bray_1981)

df_bc_all$Region <- factor(df_bc_all$Region, levels=c('NFLS', 'SS', 'NEUS'))

gg_bc_all <-  ggplot(data = df_bc_all, aes(x = Year, y = BC_simm, color = Region)) +
  geom_line(linewidth = 2) +
  #geom_point(fill = "black", size = 4, shape = 21, stroke = 0.5) +
  theme_bw(base_size = 14) +
  ylim(0.6, 1) +
  facet_wrap(~ Region) +
  scale_color_viridis_d() +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        text = element_text(family = "Arial"), 
        legend.position = "none")  

gg_bc_all

ggsave("../Scale_Dependence_Structural_Stability_LMEs/Figures/Figure SX - Bray-Curtis Similarities.jpeg", plot = gg_bc_all, width = 10, height = 3, dpi = 300)

