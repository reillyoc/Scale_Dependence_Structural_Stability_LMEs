# Dirichlet Models

# Author(s): Reilly O'Connor
# Version: 2025-08-08

# Load Pkgs
library(tidyverse)
library(vegan)
library(zoo)
library(reshape2)
library(cowplot)
library(beepr)
library(easystats)
library(bbmle)
library(DirichletReg)
library(forecast)
library(caret)

# load data
df_sp_tips <- read.csv("../Scale_Dependence_Structural_Stability_LMEs/Data/Output Data/richness_tips_split_chao.csv", header = T) %>%
  filter(year < 2004)

df_gf_biomass <- read.csv("../Scale_Dependence_Structural_Stability_LMEs/Data/Output Data/groundfish sum biomass.csv", header = T) %>%
  rename(year = Year)

df_nao <- read.csv("../Scale_Dependence_Structural_Stability_LMEs/Data/NAO Data/NAO 1973-2004.csv", header = T) %>%
  rename(year = Year) %>%
  filter(year >= 1973) 

df_temp_harvest <- read.csv("../Scale_Dependence_Structural_Stability_LMEs/Data/Cumulative Harvest Pressure Data.csv", header = T) %>%
  rename(year = Year) %>%
  dplyr::select(year, Region, total_landings_km2) %>%
  mutate(log_landings = log10(total_landings_km2)) %>%
  group_by(Region) %>%
  arrange(Region, year) %>%
  mutate(lag_landings = lag(total_landings_km2),
         lag_log_landings = lag(log_landings)) %>%
  ungroup() %>%
  filter(year >= 1973 & year < 2004)

df_nfls_harv <- df_temp_harvest %>%
  filter(Region == "NFLS")

df_ss_harv <- df_temp_harvest %>%
  filter(Region == "SS")

df_neus_harv <- df_temp_harvest %>%
  filter(Region == "NEUS")


###### Split landings into long-term trends + anomalies #####
#select best span for smoothing using k-fold cross validation
ctrl <- trainControl(method = "cv", number = 5)
grid <- expand.grid(span = seq(0.5, 0.9, len = 5), degree = 1)

#perform cross-validation using smoothing spans ranging from 0.5 to 0.9
nfls_loess_model <- train(lag_landings ~ year, data = df_nfls_harv, method = "gamLoess", tuneGrid=grid, trControl = ctrl)
nfls_loess_model

ss_loess_model <- train(lag_landings ~ year, data = df_ss_harv, method = "gamLoess", tuneGrid=grid, trControl = ctrl)
ss_loess_model

neus_loess_model <- train(lag_landings ~ year, data = df_neus_harv, method = "gamLoess", tuneGrid=grid, trControl = ctrl)
neus_loess_model

loess_nfls <- loess(lag_landings ~ year, data = df_nfls_harv, span = 0.5)
summary(loess_nfls)
df_nfls_harv$landings_trend <- predict(loess_nfls)
df_nfls_harv$landings_anom <- residuals(loess_nfls)

acf(loess_nfls$residuals)

plot(df_nfls_harv$lag_landings)
plot(df_nfls_harv$landings_trend)
plot(df_nfls_harv$landings_anom)

#Scotian Shelf
loess_ss <- loess(lag_landings ~ year, data = df_ss_harv, span = 0.5)
summary(loess_ss)

df_ss_harv$landings_trend <- predict(loess_ss)
df_ss_harv$landings_anom <- residuals(loess_ss)

acf(loess_ss$residuals)

plot(df_ss_harv$lag_landings)
plot(df_ss_harv$landings_trend)
plot(df_ss_harv$landings_anom)


#Northeast US
loess_neus <- loess(lag_landings ~ year, data = df_neus_harv, span = 0.5)
summary(loess_neus)
df_neus_harv$landings_trend <- predict(loess_neus)
df_neus_harv$landings_anom  <- residuals(loess_neus)

acf(loess_neus$residuals)

plot(df_neus_harv$lag_landings)
plot(df_neus_harv$landings_trend)
plot(df_neus_harv$landings_anom)

df_temp_harvest_trend_anom <- rbind(df_nfls_harv, df_ss_harv, df_neus_harv)

df_temp_harvest_trend_anom_diff <- df_temp_harvest_trend_anom %>%
  arrange(Region, year) %>%
  group_by(Region) %>%
  mutate(lag_landings_diff = residuals(arima(lag_landings, order = c(0,1,0)))) %>%
  ungroup()


##### Organize data, calculate scaled predictor variables, use chao corrected spp rich estimates #####
df_nfls_sp_sr <- df_sp_tips %>%
  dplyr::select(year, Trophic_interval, Region, total_sr = total_sr_chao) %>%
  filter(Trophic_interval == "3.0-3.4") %>%
  filter(Region == "NFLS") %>%
  dplyr::select(-Trophic_interval) %>%
  ungroup()

df_sp_sr <- df_sp_tips %>%
  dplyr::select(year, Trophic_interval, Region, total_sr = total_sr_chao) %>%
  filter(Trophic_interval == "3.0-3.4") %>%
  filter(! (Region == "NFLS")) %>%
  dplyr::select(-Trophic_interval) %>%
  group_by(Region) %>%
  arrange(Region, year)

df_sp_tips$group_id <- interaction(df_sp_tips$Region, df_sp_tips$Trophic_interval)

df_nfls_tips <- df_sp_tips %>%
  filter(Region == "NFLS") %>%
  dplyr::select(year, Trophic_interval, Region, tips = tips_chao)

df_sp_tips_sr <- df_sp_tips %>%
  filter(! (Region == "NFLS")) %>%
  dplyr::select(year, Trophic_interval, Region, tips = tips_chao) %>%
  rbind(df_nfls_tips) %>%
  left_join(df_nao, by = join_by(year)) %>%
  left_join(df_temp_harvest_trend_anom_diff, by = join_by(year, Region)) %>%
  left_join(df_sp_sr, by = join_by(year, Region)) %>%
  left_join(df_gf_biomass, by = join_by(year, Region))

df_sp_tips_sub <- df_sp_tips_sr %>%
  mutate(log_total_sr = log10(total_sr)) %>%
  group_by(Region) %>%
  mutate(scale_sr = scale(total_sr)[,1],
         scale_nao = scale(DJFM_NAO_Index)[,1],
         scale_anom_landings = scale(landings_anom)[,1],
         scale_landings = scale(landings_trend)[,1],
         scale_landings_diff = scale(lag_landings_diff)[,1],
         scale_year = scale(year)[,1]) %>%
  filter(! (year == 2004)) %>%
  dplyr::select(year, scale_year, Region, Trophic_Interval = Trophic_interval, tips, scale_sr, scale_nao,  scale_anom_landings, scale_landings, scale_landings_diff, scale_sum_geom_dens, DJFM_NAO_Index, landings_anom, landings_trend, lag_landings_diff) %>%
  mutate(Region = factor(Region))

df_sp_tis <- df_sp_tips_sub %>%
  pivot_wider(id_cols = c(year, Region), 
              values_from = tips, 
              names_from = Trophic_Interval)

df_sp_tis_covar <- df_sp_tips_sub %>%
  distinct(year, scale_year, Region, scale_sr, scale_nao, scale_anom_landings, scale_landings, scale_landings_diff, scale_sum_geom_dens, DJFM_NAO_Index, landings_anom, landings_trend, lag_landings_diff) %>%
  left_join(df_sp_tis, by = join_by(year, Region))

df_sp_tips_sub$Region <- factor(df_sp_tips_sub$Region, levels=c('NFLS', 'SS', 'NEUS'))

gg_tips <- ggplot(df_sp_tips_sub, aes(x = year, y = tips, color = Trophic_Interval)) +
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

# ggsave("../Scale_Dependence_Structural_Stability_LMEs/Figures/Figure 3 - Temporal TIPs.jpeg", plot = gg_tips, dpi = 300, width = 10, height = 4)

#Community NMDS & Dissimilarity Plots
source("../Scale_Dependence_Structural_Stability_LMEs/src/4 - Temporal Community Dissimilarity NMDS.R")

gg_tips_bc <- plot_grid(gg_tips, gg_bc_all, nrow = 2, align = "hv")

gg_tips_bc

# ggsave("../Scale_Dependence_Structural_Stability_LMEs/Figures/Figure 3 - Temporal TIPs + BC.jpeg", plot = gg_tips_bc, dpi = 300, width = 10, height = 8)


#Look for change over time for regions and trophic intervals
df_sp_tis_covar_nfls <- df_sp_tis_covar %>%
  filter(Region == "NFLS")

hist(df_sp_tis_covar_nfls$`3.0-3.4`)
hist(df_sp_tis_covar_nfls$`3.5-3.9`)
hist(df_sp_tis_covar_nfls$`4.0-4.4`)
hist(df_sp_tis_covar_nfls$`4.5-4.9`)

summary(lm(`3.0-3.4` ~ year, data = df_sp_tis_covar_nfls))
summary(lm(`3.5-3.9` ~ year, data = df_sp_tis_covar_nfls))
summary(lm(`4.0-4.4` ~ year, data = df_sp_tis_covar_nfls))
summary(lm(`4.5-4.9` ~ year, data = df_sp_tis_covar_nfls))

df_sp_tis_covar_ss <- df_sp_tis_covar %>%
  filter(Region == "SS")

summary(lm(`3.0-3.4` ~ year, data = df_sp_tis_covar_ss))
summary(lm(`3.5-3.9` ~ year, data = df_sp_tis_covar_ss))
summary(lm(`4.0-4.4` ~ year, data = df_sp_tis_covar_ss))
summary(lm(`4.5-4.9` ~ year, data = df_sp_tis_covar_ss))

df_sp_tis_covar_neus <- df_sp_tis_covar %>%
  filter(Region == "NEUS")


summary(lm(`3.0-3.4` ~ year, data = df_sp_tis_covar_neus))
summary(lm(`3.5-3.9` ~ year, data = df_sp_tis_covar_neus))
summary(lm(`4.0-4.4` ~ year, data = df_sp_tis_covar_neus))
summary(lm(`4.5-4.9` ~ year, data = df_sp_tis_covar_neus))


##### Start Dirichlet Model Selection Process #####
#Set SS (Middle Region) as baseline
df_sp_tis_covar$Region <- relevel(df_sp_tis_covar$Region, ref = "SS")

#Take a look at dispersion and means
df_sp_tis_covar_reg <- df_sp_tis_covar %>%
  filter(Region == "NFLS")

plot(DR_data(df_sp_tis_covar_reg[, c("3.0-3.4", "3.5-3.9", "4.0-4.4")]))

#Convert proportions to Dirichlet format
df_sp_tis_covar$Y <- DR_data(df_sp_tis_covar[, c("4.5-4.9", "4.0-4.4", "3.5-3.9", "3.0-3.4")])

##### Model selection starting with full mean model and precision models #####
#Precision model first, test interactions
mod0_all <- DirichReg(Y ~ (scale_year + scale_landings + scale_anom_landings + scale_nao)*Region |
                        (scale_year + scale_landings + scale_anom_landings + scale_nao)*Region,
                      data = df_sp_tis_covar,
                      model = "alternative", control = list(iterlim = 6000, tol1 = 1e-5, tol2 = 1e-10),
                      verbosity = 1)
summary(mod0_all)

mod0_all.1 <- DirichReg(Y ~  (scale_year + scale_landings + scale_anom_landings + scale_nao)*Region |
                          scale_landings + (scale_year + scale_anom_landings + scale_nao)*Region,
                        data = df_sp_tis_covar,
                        model = "alternative", control = list(iterlim = 6000, tol1 = 1e-6, tol2 = 1e-12),
                        verbosity = 1)
summary(mod0_all.1)

mod0_all.2 <- DirichReg(Y ~ (scale_year + scale_landings + scale_anom_landings + scale_nao)*Region |
                          scale_anom_landings + (scale_landings + scale_year + scale_nao)*Region,
                        data = df_sp_tis_covar,
                        model = "alternative", control = list(iterlim = 6000, tol1 = 1e-6, tol2 = 1e-12),
                        verbosity = 1)
summary(mod0_all.2)

mod0_all.3 <- DirichReg(Y ~ (scale_year + scale_landings + scale_anom_landings + scale_nao)*Region |
                          scale_nao + (scale_landings + scale_anom_landings + scale_year)*Region,
                        data = df_sp_tis_covar,
                        model = "alternative", control = list(iterlim = 6000, tol1 = 1e-6, tol2 = 1e-12),
                        verbosity = 1)
summary(mod0_all.3)

mod0_all.4 <- DirichReg(Y ~ (scale_year + scale_landings + scale_anom_landings + scale_nao)*Region |
                          scale_year + (scale_landings + scale_anom_landings + scale_nao)*Region,
                        data = df_sp_tis_covar,
                        model = "alternative", control = list(iterlim = 6000, tol1 = 1e-6, tol2 = 1e-12),
                        verbosity = 1)
summary(mod0_all.4)

anova.DirichletRegModel(mod0_all, mod0_all.1, mod0_all.2, mod0_all.3, mod0_all.4)

#Model fit declines only with removal of nao from region interaction...
#Remove all but NAO from region interaction
mod1_all <- DirichReg(Y ~  (scale_year + scale_landings + scale_anom_landings + scale_nao)*Region |
                        scale_landings + scale_anom_landings + scale_year + (scale_nao)*Region,
                        data = df_sp_tis_covar,
                        model = "alternative", control = list(iterlim = 6000, tol1 = 1e-5, tol2 = 1e-12),
                        verbosity = 1)
summary(mod1_all)

anova.DirichletRegModel(mod0_all, mod1_all)
AICctab(mod0_all, mod1_all)

#mod1_all similar fit but AICc much improved 
mod1_all.1 <- DirichReg(Y ~ scale_year + (scale_landings + scale_anom_landings + scale_nao)*Region |
                          scale_landings + scale_anom_landings + scale_year + (scale_nao)*Region,
                      data = df_sp_tis_covar,
                      model = "alternative", control = list(iterlim = 6000, tol1 = 1e-6, tol2 = 1e-12),
                      verbosity = 1)
summary(mod1_all.1)

mod1_all.2 <- DirichReg(Y ~  scale_landings + (scale_year + scale_anom_landings + scale_nao)*Region |
                          scale_landings + scale_anom_landings + scale_year + (scale_nao)*Region,
                        data = df_sp_tis_covar,
                        model = "alternative", control = list(iterlim = 6000, tol1 = 1e-6, tol2 = 1e-12),
                        verbosity = 1)
summary(mod1_all.2)

mod1_all.3 <- DirichReg(Y ~  scale_anom_landings + (scale_landings + scale_year + scale_nao)*Region |
                          scale_landings + scale_anom_landings + scale_year + (scale_nao)*Region,
                        data = df_sp_tis_covar,
                        model = "alternative", control = list(iterlim = 6000, tol1 = 1e-5, tol2 = 1e-12),
                        verbosity = 1)
summary(mod1_all.3)

mod1_all.4 <- DirichReg(Y ~  scale_nao + (scale_landings + scale_year + scale_anom_landings)*Region |
                          scale_landings + scale_anom_landings + scale_year + (scale_nao)*Region,
                        data = df_sp_tis_covar,
                        model = "alternative", control = list(iterlim = 6000, tol1 = 1e-6, tol2 = 1e-12),
                        verbosity = 1)
summary(mod1_all.4)

anova.DirichletRegModel(mod1_all, mod1_all.1, mod1_all.2, mod1_all.3, mod1_all.4)

#removal of year with region interaction nearly significant... keep in interaction
mod1_all.5 <- DirichReg(Y ~  scale_nao + scale_landings + scale_anom_landings + (scale_year)*Region |
                          scale_landings + scale_anom_landings + scale_year + (scale_nao)*Region,
                        data = df_sp_tis_covar,
                        model = "alternative", control = list(iterlim = 6000, tol1 = 1e-6, tol2 = 1e-12),
                        verbosity = 1)
summary(mod1_all.5)

anova.DirichletRegModel(mod1_all, mod1_all.1, mod1_all.2, mod1_all.3, mod1_all.4, mod1_all.5)
AICctab(mod1_all, mod1_all.1, mod1_all.2, mod1_all.3, mod1_all.4, mod1_all.5)

#Drop interactions in mean model except for yearXregion
#Test whether to drop variables in precision equation
mod2_all <-  DirichReg(Y ~  scale_nao + scale_landings + scale_anom_landings + (scale_year)*Region |
                         scale_landings + scale_anom_landings + scale_year + (scale_nao)*Region,
                       data = df_sp_tis_covar,
                       model = "alternative", control = list(iterlim = 6000, tol1 = 1e-6, tol2 = 1e-12),
                       verbosity = 1)
summary(mod2_all)

mod2_all.1 <-  DirichReg(Y ~  scale_nao + scale_landings + scale_anom_landings + (scale_year)*Region |
                           scale_landings + scale_anom_landings + (scale_nao)*Region,
                         data = df_sp_tis_covar,
                         model = "alternative", control = list(iterlim = 6000, tol1 = 1e-6, tol2 = 1e-10),
                         verbosity = 1)
summary(mod2_all.1)

mod2_all.2 <-  DirichReg(Y ~  scale_nao + scale_landings + scale_anom_landings + (scale_year)*Region |
                           scale_landings + scale_year + (scale_nao)*Region,
                         data = df_sp_tis_covar,
                         model = "alternative", control = list(iterlim = 6000, tol1 = 1e-6, tol2 = 1e-10),
                         verbosity = 1)
summary(mod2_all.2)

mod2_all.3 <-  DirichReg(Y ~  scale_nao + scale_landings + scale_anom_landings + (scale_year)*Region |
                           scale_year + scale_anom_landings + (scale_nao)*Region,
                         data = df_sp_tis_covar,
                         model = "alternative", control = list(iterlim = 6000, tol1 = 1e-6, tol2 = 1e-10),
                         verbosity = 1)
summary(mod2_all.3)

anova.DirichletRegModel(mod2_all, mod2_all.1, mod2_all.2, mod2_all.3)
AICctab(mod2_all, mod2_all.1, mod2_all.2, mod2_all.3)

#Drop variables from mean equation, landings and year model similar fits
#Select landings model as it is a measured variable and a slightly bettter fit than year
#Drop year...
mod3_all <-  DirichReg(Y ~ scale_nao + scale_landings + scale_anom_landings + (scale_year)*Region |
                         scale_landings + scale_anom_landings + (scale_nao)*Region,
                         data = df_sp_tis_covar,
                         model = "alternative", control = list(iterlim = 6000, tol1 = 1e-6, tol2 = 1e-10),
                         verbosity = 1)
summary(mod3_all)

mod3_all.1 <-  DirichReg(Y ~ scale_nao + scale_landings + (scale_year)*Region |
                         scale_landings + scale_anom_landings + (scale_nao)*Region,
                       data = df_sp_tis_covar,
                       model = "alternative", control = list(iterlim = 6000, tol1 = 1e-6, tol2 = 1e-10),
                       verbosity = 1)
summary(mod3_all.1)

mod3_all.2 <-  DirichReg(Y ~ scale_nao + scale_anom_landings + (scale_year)*Region |
                         scale_landings + scale_anom_landings + (scale_nao)*Region,
                       data = df_sp_tis_covar,
                       model = "alternative", control = list(iterlim = 6000, tol1 = 1e-6, tol2 = 1e-10),
                       verbosity = 1)
summary(mod3_all.2)

mod3_all.3 <-  DirichReg(Y ~ scale_landings + scale_anom_landings + (scale_year)*Region |
                         scale_landings + scale_anom_landings + (scale_nao)*Region,
                       data = df_sp_tis_covar,
                       model = "alternative", control = list(iterlim = 6000, tol1 = 1e-6, tol2 = 1e-10),
                       verbosity = 1)
summary(mod3_all.3)

anova.DirichletRegModel(mod3_all, mod3_all.1, mod3_all.2, mod3_all.3)
AICctab(mod3_all, mod3_all.1, mod3_all.2, mod3_all.3)

summary(mod3_all.1)

#landing and landings anomalies no significant effect on means, remove and re-evaulate mean + precision models without them in mean
#Finish optimizing mean model
mod4_all <-  DirichReg(Y ~ scale_nao + scale_year*Region |
                         scale_landings + scale_anom_landings + (scale_nao)*Region,
                       data = df_sp_tis_covar,
                       model = "alternative", control = list(iterlim = 6000, tol1 = 1e-6, tol2 = 1e-10),
                       verbosity = 1)
summary(mod4_all)

mod4_all.1 <-  DirichReg(Y ~ scale_nao + scale_year + Region |
                         scale_landings + scale_anom_landings + (scale_nao)*Region,
                       data = df_sp_tis_covar,
                       model = "alternative", control = list(iterlim = 6000, tol1 = 1e-6, tol2 = 1e-10),
                       verbosity = 1)
summary(mod4_all.1)

mod4_all.2 <-  DirichReg(Y ~ scale_nao + Region |
                           scale_landings + scale_anom_landings + (scale_nao)*Region,
                         data = df_sp_tis_covar,
                         model = "alternative", control = list(iterlim = 6000, tol1 = 1e-6, tol2 = 1e-10),
                         verbosity = 1)
summary(mod4_all.2)

mod4_all.3 <-  DirichReg(Y ~ scale_year + Region |
                           scale_landings + scale_anom_landings + (scale_nao)*Region,
                         data = df_sp_tis_covar,
                         model = "alternative", control = list(iterlim = 6000, tol1 = 1e-6, tol2 = 1e-10),
                         verbosity = 1)
summary(mod4_all.3)

anova.DirichletRegModel(mod3_all.2, mod4_all, mod4_all.1, mod4_all.2, mod4_all.3)
AICctab(mod3_all.2, mod4_all, mod4_all.1, mod4_all.2, mod4_all.3)

#Finish by optimizing precision model
mod4_all.4 <-  DirichReg(Y ~ scale_nao + Region |
                           scale_landings + (scale_nao)*Region,
                         data = df_sp_tis_covar,
                         model = "alternative", control = list(iterlim = 6000, tol1 = 1e-6, tol2 = 1e-10),
                         verbosity = 1)
summary(mod4_all.4)

mod4_all.5 <-  DirichReg(Y ~ scale_nao + Region |
                           scale_anom_landings + (scale_nao)*Region,
                         data = df_sp_tis_covar,
                         model = "alternative", control = list(iterlim = 6000, tol1 = 1e-6, tol2 = 1e-10),
                         verbosity = 1)
summary(mod4_all.5)

mod4_all.6 <-  DirichReg(Y ~ scale_nao + Region |
                           scale_landings + scale_anom_landings + scale_nao + Region,
                         data = df_sp_tis_covar,
                         model = "alternative", control = list(iterlim = 6000, tol1 = 1e-6, tol2 = 1e-10),
                         verbosity = 1)
summary(mod4_all.6)

anova.DirichletRegModel(mod4_all.2, mod4_all.4, mod4_all.5, mod4_all.6)

mod4_all.7 <-  DirichReg(Y ~ Region |
                           scale_landings + scale_anom_landings + scale_nao + Region,
                         data = df_sp_tis_covar,
                         model = "alternative", control = list(iterlim = 6000, tol1 = 1e-6, tol2 = 1e-10),
                         verbosity = 1)
summary(mod4_all.7)

anova.DirichletRegModel(mod4_all.2, mod4_all.7)

#mod4_all.2 is the best fit model after model selection
#mod4_all.4 is similar fit, but given anova suggests that removing scale_anom_landings significantly reduces fit
#Results are similar between models
##### Plot Final Model #####
mod_final <- mod4_all.2
summary.DirichletRegModel(mod_final)
confint.DirichletRegModel(mod_final, level = 0.90)

summ_mod3_all <- summary.DirichletRegModel(mod_final)

mod_final_resid <- residuals(mod_final, type = "standardized")
df_mod_final_resid <- as.data.frame(mod_final_resid)

#look at residuals
plot(df_mod_final_resid$x[,"4.5-4.9"])
plot(df_mod_final_resid$x[,"4.0-4.4"])
plot(df_mod_final_resid$x[,"3.5-3.9"])
plot(df_mod_final_resid$x[,"3.0-3.4"])

#look at residual acfs
acf(df_mod_final_resid$x[,"4.5-4.9"])
acf(df_mod_final_resid$x[,"4.0-4.4"])
acf(df_mod_final_resid$x[,"3.5-3.9"])
acf(df_mod_final_resid$x[,"3.0-3.4"])


##### Plot Mean Effects #####
# Predicted mean proportions
pred_mean <- predict(mod_final, type = "response")

df_pred_mean <- data.frame(
  year = df_sp_tis_covar$year,
  Region = df_sp_tis_covar$Region,
  mean_prop = pred_mean
)

df_long_mean <- df_pred_mean %>%
  # turn the four mean_prop.* columns into two columns: TIS + mean_prop
  pivot_longer(
    cols = starts_with("mean_prop."),
    names_to  = "TIS_raw",
    values_to = "mean_prop"
  ) %>%
  mutate(
    Trophic_Interval = case_when(
      TIS_raw == "mean_prop.3.0.3.4" ~ "3.0–3.4",
      TIS_raw == "mean_prop.3.5.3.9" ~ "3.5–3.9",
      TIS_raw == "mean_prop.4.0.4.4" ~ "4.0–4.4",
      TIS_raw == "mean_prop.4.5.4.9" ~ "4.5–4.9",
      TRUE ~ NA_character_  # optional for unmatched values
    ),
    Trophic_Interval = factor(Trophic_Interval, levels = c("3.0–3.4","3.5–3.9","4.0–4.4","4.5–4.9"))
  ) %>%
  dplyr::select(year, Region, Trophic_Interval, mean_prop)

df_sp_tis_covar_long <- df_sp_tis_covar %>%
  dplyr::select(year, Region, `3.0-3.4`:`4.5-4.9`) %>%
  pivot_longer(
  cols = c(`3.0-3.4`:`4.5-4.9`),
  names_to  = "Trophic_Interval",
  values_to = "obs_prop"
)

#Regional Grounfish Biomass Collapses
df_gfb_windows <- tribble(
  ~Region, ~xintercept,
  "NFLS",  1990,
  "NFLS",  1993,
  "SS",    1992,
  "SS",    1998,
  "NEUS",  1981,
  "NEUS",  1991,
)

df_long_mean$Region <- factor(df_long_mean$Region, levels=c('NFLS', 'SS', 'NEUS'))
df_gfb_windows$Region <- factor(df_gfb_windows$Region, levels=c('NFLS', 'SS', 'NEUS'))

# Plot predicted mean proportions
gg_pred_tips <- ggplot(df_long_mean, aes(x = year, y = mean_prop, color = Trophic_Interval)) +
  geom_vline(data = df_gfb_windows, aes(xintercept = xintercept), linetype = 'dashed', alpha = 0.8) +
  geom_line(linewidth = 2, alpha = 0.8) +
  facet_wrap(~ Region) +
  ylab("Predicted mean proportion") +
  theme_bw(base_size = 14) +
  xlim(1970, 2005) +
  ylim(0, 0.5) +
  scale_fill_manual(values = c("#4245c4", "#f78c2f",  "#649f4d", "#b61790")) +
  scale_color_manual(values = c("#4245c4", "#f78c2f",  "#649f4d", "#b61790")) +
  labs(color = "Trophic Interval")

gg_pred_tips

#nao effect on Mean TIPs
pred_data_nao <- expand.grid(
  scale_nao = seq(-2, 2, length.out = 100),
  scale_year = mean(df_sp_tis_covar$scale_year),
  scale_landings    = mean(df_sp_tis_covar$scale_landings),
  scale_anom_landings    = mean(df_sp_tis_covar$scale_anom_landings),
  Region         = levels(df_sp_tis_covar$Region)
)

# φ on natural scale
mu_hat_nao <- data.frame(predict(
  mod_final,
  newdata = pred_data_nao,
  mu = TRUE, alpha = FALSE, phi = FALSE
))

pred_data_nao_mu <- pred_data_nao %>%
  cbind(mu_hat_nao) %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to  = "TIs",
    values_to = "mean_prop"
  ) %>%
  mutate(
  Trophic_Interval = case_when(
    TIs == "X1" ~ "3.0–3.4",
    TIs == "X2" ~ "3.5–3.9",
    TIs == "X3" ~ "4.0–4.4",
    TIs == "X4" ~ "4.5–4.9",
    TRUE ~ NA_character_  # optional for unmatched values
  ),
  Trophic_Interval = factor(Trophic_Interval, levels = c("3.0–3.4","3.5–3.9","4.0–4.4","4.5–4.9"))) %>%
  dplyr::select(scale_nao, Region, Trophic_Interval, scale_landings, mean_prop)
  
pred_data_nao_mu$Region <- factor(pred_data_nao_mu$Region, levels=c('NFLS', 'SS', 'NEUS'))

gg_mu_nao <- ggplot(pred_data_nao_mu %>% filter(Region == "NFLS"), aes(scale_nao, mean_prop, color = Trophic_Interval)) +
  geom_line(linewidth = 2, alpha = 0.8) +
  # facet_wrap(~ Region) +
  ylim(0, 0.5) +
  labs(x = "NAO (Scaled)", y = "Mean Trophic Interval Proportion") +
  scale_color_manual(values = c( "#b61790", "#649f4d", "#f78c2f", "#4245c4")) +
  theme_bw(14) +
  theme(legend.position = "none")

gg_mu_nao

##### Look at Precision #####
#Predicted means (fitted values)
pred_df <- as.data.frame(predict(mod_final, se.fit = TRUE))

#Phi predicted values
pred.phi <- predict(mod_final, phi = TRUE)$phi
pred.phi <- data.frame(df_sp_tis_covar[, c(1:3)], phi = pred.phi)

df_pred.phi <- pred.phi %>%
  mutate(variable = "precision",
         variance = 1 / (phi + 1),
         log_phi = log10(phi),
         inv_phi = 1/phi,
         inv_log_phi = 1/log_phi,
         Region = as.factor(Region)) %>%
  dplyr::select(year, Region, variable, phi:inv_log_phi) %>%
  ungroup()
  

#Regional Grounfish Biomass Collapses
df_gfb_windows <- tribble(
  ~Region, ~xintercept,
  "NFLS",  1990,
  "NFLS",  1993,
  "SS",    1992,
  "SS",    1998,
  "NEUS",  1981,
  "NEUS",  1991,
)

df_pred.phi$Region <- factor(df_pred.phi$Region, levels=c('NFLS', 'SS', 'NEUS'))
df_gfb_windows$Region <- factor(df_gfb_windows$Region, levels=c('NFLS', 'SS', 'NEUS'))

gg_pred_precis <- ggplot(df_pred.phi, aes(x = year, y = phi, color = Region)) +
  geom_vline(data = df_gfb_windows, aes(xintercept = xintercept), linetype = "dashed", alpha = 0.8) +
  geom_line(linewidth = 2, alpha = 0.8) +
  facet_wrap(~Region, nrow = 1, scale = "free_y") +
  theme_bw(base_size = 14) +
  scale_color_viridis_d() +
  scale_y_log10() +
  xlim(1970, 2006) +
  labs(color = "Region",
       y = expression(phi),
       x = "Year") +
  scale_y_continuous(expand = expansion(mult = c(0.25, 0.25)))

gg_pred_precis

gg_fig4 <- cowplot::plot_grid(gg_pred_tips, gg_pred_precis, ncol = 1, align = "hv")

gg_fig4

ggsave("../Scale_Dependence_Structural_Stability_LMEs/Figures/Figure 4 - Model Estimates - Mean + Precision.jpeg", plot = gg_fig4, width = 10, height = 6, dpi = 300)

ggsave("../Scale_Dependence_Structural_Stability_LMEs/Figures/Figure 4 - Model Estimates - Precision.jpeg", plot = gg_pred_precis, width = 10, height = 3, dpi = 300)

ggsave("../Scale_Dependence_Structural_Stability_LMEs/Figures/Figure 4 - Model Estimates - Mean.jpeg", plot = gg_pred_tips, width = 10, height = 3, dpi = 300)

#Landings Effect on Precision (Variability in Trophic Structure)
pred_data_land <- expand.grid(
  scale_anom_landings = seq(-2, 2, length.out = 100),
  # scale_year     = mean(df_sp_tis_covar$scale_year),
  scale_landings = mean(df_sp_tis_covar$scale_landings),
  scale_nao = mean(df_sp_tis_covar$scale_nao),
  Region         = unique(df_sp_tis_covar$Region)   # keep factor levels
)

# φ on natural scale
phi_hat_land <- as.numeric(predict(
  mod_final,
  newdata = pred_data_land,
  mu = FALSE, alpha = FALSE, phi = TRUE
))

# log(φ) = link scale
pred_data_land$phi     <- phi_hat_land
pred_data_land$log_phi <- log10(phi_hat_land)

pred_data_land$Region <- factor(pred_data_land$Region, levels=c('NFLS', 'SS', 'NEUS'))

gg_phi_land <- ggplot(pred_data_land, aes(scale_anom_landings, phi, color = Region)) +
  geom_line(linewidth = 2, alpha = 0.8) +
  # facet_wrap(~ Region, scale = "free") +
  labs(x = "Scaled landings", y = expression(1/phi)) +
  scale_color_viridis_d() +
  theme_bw(14) +
  scale_y_log10(limits = c(25, 600)) +
  theme(legend.position = "none")

gg_phi_land

#Landings Anomalies Effect on Precision (Variability in Trophic Structure)
pred_data_land_anom <- expand.grid(
  scale_anom_landings = seq(-2, 2, length.out = 100),
  scale_landings     = mean(df_sp_tis_covar$scale_landings),
  scale_nao = mean(df_sp_tis_covar$scale_nao),
  Region         = unique(df_sp_tis_covar$Region)   # keep factor levels
)

# φ on natural scale
phi_hat_land_anom <- as.numeric(predict(
  mod_final,
  newdata = pred_data_land_anom,
  mu = FALSE, alpha = FALSE, phi = TRUE
))

# log(φ) = link scale
pred_data_land_anom$phi     <- phi_hat_land_anom
pred_data_land_anom$log_phi <- log10(phi_hat_land_anom)

pred_data_land_anom$Region <- factor(pred_data_land_anom$Region, levels=c('NFLS', 'SS', 'NEUS'))

gg_phi_land_anom <- ggplot(pred_data_land_anom, aes(scale_anom_landings, phi, color = Region)) +
  geom_line(linewidth = 2, alpha = 0.8) +
  # facet_wrap(~ Region, scale = "free") +
  labs(x = "Scaled Landings Anomalies", y = expression(log(phi))) +
  scale_color_viridis_d() +
  theme_bw(14) +
  scale_y_log10(limits = c(25, 600)) +
  theme(legend.position = "none",
        axis.title.x = element_blank())

gg_phi_land_anom

#NAO Precision (Variability in Trophic Structure)
pred_data_nao <- expand.grid(
  scale_nao = seq(-2, 2, length.out = 100),
  scale_landings    = mean(df_sp_tis_covar$scale_landings),
  scale_anom_landings    = mean(df_sp_tis_covar$scale_anom_landings),
  Region         = levels(df_sp_tis_covar$Region)   # keep factor levels
)

# φ on natural scale
phi_hat_nao <- as.numeric(predict(
  mod_final,
  newdata = pred_data_nao,
  mu = FALSE, alpha = FALSE, phi = TRUE
))

# log(φ) = link scale
pred_data_nao$phi     <- phi_hat_nao
pred_data_nao$log_phi <- log10(phi_hat_nao)

pred_data_nao$Region <- factor(pred_data_nao$Region, levels=c('NFLS', 'SS', 'NEUS'))

gg_phi_nao <- ggplot(pred_data_nao, aes(scale_nao, phi, color = Region)) +
  geom_line(linewidth = 2, alpha = 0.8) +
  facet_wrap(~ Region, ncol = 1) +
  labs(x = "North Atlantic Oscillation (Scaled)", y = expression(1/phi)) +
  scale_color_viridis_d() +
  # ylim(0.0, 0.025) +
  scale_y_log10(limits = c(25, 600)) +
  theme_bw(14) +
  theme(legend.position = "none",
        axis.title.x = element_blank())

gg_phi_nao

gg_phi_land_anoms <- cowplot::plot_grid(gg_phi_land, gg_phi_land_anom, nrow = 2, align = "hv")
gg_phi_land_anoms

gg_phi_land_nao <- cowplot::plot_grid(gg_phi_land_anoms, gg_phi_nao, nrow = 1, align = "hv")

gg_phi_land_nao

ggsave("../Scale_Dependence_Structural_Stability_LMEs/Figures/Figure 5 - TIPs Phi Panel Landings NAO.jpeg", plot = gg_phi_land_nao, width = 8, height = 5, dpi = 300)

