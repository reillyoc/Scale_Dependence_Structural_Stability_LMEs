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

lm_nfls <- lm(lag_landings ~ year, data = df_nfls_harv)
summary(lm_nfls)

df_nfls_harv$landings_trend <- predict(loess_nfls)
df_nfls_harv$landings_anom <- residuals(loess_nfls)
df_nfls_harv$landings_detrend <- residuals(lm_nfls)

acf(loess_nfls$residuals)
acf(lm_nfls$residuals)

plot(df_nfls_harv$lag_landings)
plot(df_nfls_harv$landings_trend)
plot(df_nfls_harv$landings_anom)

#Scotian Shelf
loess_ss <- loess(lag_landings ~ year, data = df_ss_harv, span = 0.5)
summary(loess_ss)

lm_ss <- lm(lag_landings ~ year, data = df_ss_harv)
summary(lm_ss)

df_ss_harv$landings_trend <- predict(loess_ss)
df_ss_harv$landings_anom <- residuals(loess_ss)
df_ss_harv$landings_detrend <- residuals(lm_ss)

acf(loess_ss$residuals)

plot(df_ss_harv$lag_landings)
plot(df_ss_harv$landings_trend)
plot(df_ss_harv$landings_anom)


#Northeast US
loess_neus <- loess(lag_landings ~ year, data = df_neus_harv, span = 0.5)
summary(loess_neus)

lm_neus <- lm(lag_landings ~ year, data = df_neus_harv)
summary(lm_neus)

df_neus_harv$landings_trend <- predict(loess_neus)
df_neus_harv$landings_anom  <- residuals(loess_neus)
df_neus_harv$landings_detrend <- residuals(lm_neus)

acf(loess_neus$residuals)

plot(df_neus_harv$lag_landings)
plot(df_neus_harv$landings_trend)
plot(df_neus_harv$landings_anom)
plot(df_neus_harv$landings_detrend)
df_temp_harvest_trend_anom <- rbind(df_nfls_harv, df_ss_harv, df_neus_harv)

df_temp_harvest_trend_anom_diff <- df_temp_harvest_trend_anom %>%
  arrange(Region, year) %>%
  group_by(Region) %>%
  mutate(lag_landings_diff = residuals(arima(lag_landings, order = c(0,1,0)))) %>%
  ungroup()

df_temp_harvest_trend_anom_diff$Region <- factor(df_temp_harvest_trend_anom_diff$Region, levels=c('NFLS', 'SS', 'NEUS'))

gg_lag_landings <- ggplot(df_temp_harvest_trend_anom_diff, aes(x = year, y = lag_landings, color = Region)) +
  geom_line(linewidth = 2) +
  facet_wrap(~ Region, scale = "free_y") + 
  xlim(1970, 2005) +
  scale_color_viridis_d() +
  theme_bw(base_size = 14)

gg_lag_landings

gg_trend_landings <- ggplot(df_temp_harvest_trend_anom_diff, aes(x = year, y = landings_trend, color = Region)) +
  geom_line(linewidth = 2) +
  facet_wrap(~ Region, scale = "free_y") + 
  xlim(1970, 2005) +
  scale_color_viridis_d() +
  theme_bw(base_size = 14)

gg_trend_landings

gg_anom_landings <- ggplot(df_temp_harvest_trend_anom_diff, aes(x = year, y = landings_anom, color = Region)) +
  geom_line(linewidth = 2) +
  facet_wrap(~ Region, scale = "free_y") + 
  xlim(1970, 2005) +
  scale_color_viridis_d() +
  theme_bw(base_size = 14)

gg_anom_landings

gg_landings_grid <- cowplot::plot_grid(gg_lag_landings, gg_trend_landings, gg_anom_landings,
                                       nrow = 3, align = "hv")

gg_landings_grid

# ggsave("../Scale_Dependence_Structural_Stability_LMEs/Figures/Figure SX - Landings Decomposition.jpeg", plot = gg_landings_grid, dpi = 300, width = 10, height = 9)


##### Organize data, calculate scaled predictor variables, use chao corrected spp rich estimates #####
df_nfls_sp_sr <- df_sp_tips %>%
  dplyr::select(year, Trophic_interval, Region, total_sr = total_sr_jack2) %>%
  filter(Trophic_interval == "3.0-3.4") %>%
  filter(Region == "NFLS") %>%
  dplyr::select(-Trophic_interval) %>%
  ungroup()

df_sp_sr <- df_sp_tips %>%
  dplyr::select(year, Trophic_interval, Region, total_sr = total_sr_jack2) %>%
  filter(Trophic_interval == "3.0-3.4") %>%
  filter(! (Region == "NFLS")) %>%
  dplyr::select(-Trophic_interval) %>%
  group_by(Region) %>%
  arrange(Region, year) %>%
  rbind(df_nfls_sp_sr)

df_sp_tips$group_id <- interaction(df_sp_tips$Region, df_sp_tips$Trophic_interval)

df_nfls_tips <- df_sp_tips %>%
  filter(Region == "NFLS") %>%
  dplyr::select(year, Trophic_interval, Region, tips = tips_jack2)

df_sp_tips_sr <- df_sp_tips %>%
  filter(! (Region == "NFLS")) %>%
  dplyr::select(year, Trophic_interval, Region, tips = tips_jack2) %>%
  rbind(df_nfls_tips) %>%
  left_join(df_nao, by = join_by(year)) %>%
  left_join(df_temp_harvest_trend_anom_diff, by = join_by(year, Region)) %>%
  left_join(df_sp_sr, by = join_by(year, Region)) %>%
  left_join(df_gf_biomass, by = join_by(year, Region))

mean_sr <- df_sp_tips_sr %>%
  filter(Trophic_interval == '3.0-3.4') %>%
  group_by(Region) %>%
  summarise(mean_sr = mean(total_sr),
            se_sr = standard_error(total_sr))

df_sp_tips_sub <- df_sp_tips_sr %>%
  mutate(log_total_sr = log10(total_sr)) %>%
  group_by(Region) %>%
  mutate(scale_sr = scale(total_sr)[,1],
         scale_nao = scale(DJFM_NAO_Index)[,1],
         scale_landings = scale(lag_log_landings)[,1],
         scale_anom_landings = scale(landings_anom)[,1],
         scale_landings_trend = scale(landings_trend)[,1],
         scale_landings_diff = scale(lag_landings_diff)[,1],
         scale_landings_detrend = scale(landings_detrend)[,1],
         scale_year = scale(year)[,1]) %>%
  filter(! (year == 2004)) %>%
  dplyr::select(year, scale_year, Region, Trophic_Interval = Trophic_interval, tips, scale_sr, scale_nao,  scale_anom_landings, scale_landings, scale_landings_detrend, scale_landings_trend,  scale_landings_diff, scale_sum_geom_dens, DJFM_NAO_Index, landings_anom, landings_trend, lag_landings_diff, total_sr) %>%
  mutate(Region = factor(Region))

df_sp_tis <- df_sp_tips_sub %>%
  pivot_wider(id_cols = c(year, Region), 
              values_from = tips, 
              names_from = Trophic_Interval)

df_sp_tis_covar <- df_sp_tips_sub %>%
  distinct(year, scale_year, Region, scale_sr, scale_nao, scale_anom_landings, scale_landings, scale_landings_detrend, scale_landings_trend, scale_landings_diff, scale_sum_geom_dens, DJFM_NAO_Index, landings_anom, landings_trend, lag_landings_diff, total_sr) %>%
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
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  scale_x_continuous(limits = c(1970, 2005), breaks = scales::pretty_breaks()) +
  facet_wrap(~ Region) +
  theme(text = element_text(family = "Arial"))

gg_tips

# ggsave("../Scale_Dependence_Structural_Stability_LMEs/Figures/Figure 3 - Temporal TIPs.jpeg", plot = gg_tips, dpi = 300, width = 10, height = 4)

#Community NMDS & Dissimilarity Plots
source("../Scale_Dependence_Structural_Stability_LMEs/src/4 - Temporal Community Dissimilarity NMDS.R")

gg_tips_bc <- plot_grid(gg_tips, gg_bc_all, nrow = 2, align = "hv")

gg_tips_bc

# ggsave("../Scale_Dependence_Structural_Stability_LMEs/Figures/Figure 3 - Temporal TIPs + BC.jpeg", plot = gg_tips_bc, dpi = 300, width = 10, height = 8)

df_sp_tips_sub

gg_sr <- ggplot(df_sp_tips_sub %>% filter(Trophic_Interval == "3.0-3.4"), aes(x = year, y = total_sr, color = Region)) +
  geom_line(alpha = 0.8, linewidth = 2) + 
  labs(x = "Year",
       y = "Jack2 Esimated Total Species Richness") +
  theme_bw(base_size = 14) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  scale_x_continuous(limits = c(1970, 2005), breaks = scales::pretty_breaks()) +
  facet_wrap(~ Region, scale = "free_y") +
  scale_color_viridis_d() +
  theme(text = element_text(family = "Arial"))

gg_sr

gg_tips_sr <- plot_grid(gg_tips, gg_sr, nrow = 2, align = "hv")

gg_tips_sr

# ggsave("../Scale_Dependence_Structural_Stability_LMEs/Figures/Figure 3 - Temporal TIPs + SR.jpeg", plot = gg_tips_sr, dpi = 300, width = 10, height = 6)


##### Start Dirichlet Model Selection Process #####
#Set seed
set.seed(444)

#Set SS (Middle Region) as baseline
df_sp_tis_covar$Region <- relevel(df_sp_tis_covar$Region, ref = "SS")

#Take a look at dispersion and means
df_sp_tis_covar_reg <- df_sp_tis_covar %>%
  filter(Region == "NFLS")

plot(DR_data(df_sp_tis_covar_reg[, c("3.0-3.4", "3.5-3.9", "4.0-4.4")]))

#Convert proportions to Dirichlet format
df_sp_tis_covar$Y <- DR_data(df_sp_tis_covar[, c("4.5-4.9", "4.0-4.4", "3.5-3.9", "3.0-3.4")])

##### Model selection starting with full mean model and precision models #####
# Required: rerun var_subsets so they're fresh list elements
vars <- c("scale_year", 
          # "scale_landings",
          "scale_landings_trend",
          "scale_anom_landings",
          # "scale_landings_diff",
          # "scale_landings_detrend",
          "scale_nao")

get_var_subsets <- function(vars) {
  unlist(lapply(0:length(vars), function(i) combn(vars, i, simplify = FALSE)), recursive = FALSE)
}

make_hybrid_rhs <- function(var_subset, allow_region = TRUE) {
  if (length(var_subset) == 0) {
    return(if (allow_region) c("1", "Region") else "1")
  }
  
  idx <- seq_along(var_subset)
  splits <- unlist(lapply(0:length(idx), function(k) {
    lapply(combn(idx, k, simplify = FALSE), function(interact_idx) {
      list(
        main_vars     = var_subset[setdiff(idx, interact_idx)],
        interact_vars = var_subset[interact_idx]
      )
    })
  }), recursive = FALSE)
  
  out <- sapply(splits, function(sp) {
    parts <- c()
    
    # main effects
    if (length(sp$main_vars) > 0) {
      parts <- c(parts, paste(sp$main_vars, collapse = " + "))
    }
    
    # interaction effects
    if (length(sp$interact_vars) > 0) {
      parts <- c(parts, paste0("(", paste(sp$interact_vars, collapse = " + "), ")*Region"))
    }
    
    # add Region only if allow_region = TRUE AND not already implied
    if (allow_region) {
      already_in <- any(grepl("\\*Region", parts))
      if (!already_in) parts <- c(parts, "Region")
    }
    
    paste(parts, collapse = " + ")
  }, simplify = TRUE)
  
  # deduplicate across splits
  unique(out)
}

var_subsets <- get_var_subsets(vars)
rhs_combos <- expand.grid(mean_idx = seq_along(var_subsets),
                          precision_idx = seq_along(var_subsets),
                          stringsAsFactors = FALSE)

all_formulas <- do.call(c, mapply(function(m_idx, p_idx) {
  mean_vars <- var_subsets[[m_idx]]
  precision_vars <- var_subsets[[p_idx]]
  
  mean_rhss <- make_hybrid_rhs(mean_vars)
  precision_rhss <- make_hybrid_rhs(precision_vars)
  
  outer(mean_rhss, precision_rhss, Vectorize(function(m, p) {
    as.formula(paste("Y ~", m, "|", p))
  }))
}, rhs_combos$mean_idx, rhs_combos$precision_idx,
SIMPLIFY = FALSE))

# unique_formulas <- unique(all_formulas)
unique_formulas <- all_formulas[!duplicated(sapply(all_formulas, deparse))]

# Detect number of cores
num_cores <- detectCores() - 1  # leave one core free

# Parallel model fitting
models_interaction <- mclapply(unique_formulas, function(f) {
  # print(f)
  eval(bquote(DirichReg(.(f), data = df_sp_tis_covar, 
                        model = "alternative", 
                        control = list(iterlim = 6000, tol1 = 1e-5, tol2 = 1e-10))))
}, mc.cores = num_cores)

#Save models
# saveRDS(models_interaction, file = "../Scale_Dependence_Structural_Stability_LMEs/Data/Output Data/dirichlet_models.rds", compress = "xz")

# Load models
# models_interaction <- readRDS("../Scale_Dependence_Structural_Stability_LMEs/Data/Output Data/dirichlet_models.rds")

df_model_selection <- data.frame(
  Model_ID = seq_along(models_interaction),
  AIC = sapply(models_interaction, AIC),
  AICc = sapply(models_interaction, MuMIn::AICc)
)


df_model_selection_fin <- df_model_selection %>%
  mutate(AIC_diff = min(AIC) - AIC,
         AICc_diff = min(AICc) - AICc) %>%
  filter(AICc_diff >= -6) %>%
  mutate(
    AICc_delta = AICc - min(AICc),
    Weight = exp(-0.5 * AICc_delta) / sum(exp(-0.5 * AICc_delta)),
    perc_weight = Weight*100)

###### Subset Top Candidate Models based on AICc #####
df_model_selection_fin %>%
  filter(AICc_diff >= -2.6) %>%
  arrange(AICc)

##### Check Candidate Models & Summaries #####
models_interaction[[2134]]$call
summary(models_interaction[[1478]])


##### Plot Top Model #####
mod_final <- models_interaction[[2134]]
summary.DirichletRegModel(mod_final)
confint.DirichletRegModel(mod_final, level = 0.85)

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
#Simulate Data
means <- unlist(coef(mod_final))
vc <- vcov(mod_final)
N <- 10000
rnd <- rmvnorm(N, mean = means, sigma = vc)

#Generate linear predictors
modobj <- mod_final

#Design matrix (mean model)
dm <- do.call(cbind, modobj$X)  

#Design matrix (precision model)
dp <- modobj$Z                  
n <- nrow(df_sp_tis_covar)

#Extract linear predictors for each Trophic Interval
tip_3   <- rnd[, grep("3.0-3.4", colnames(rnd))]   %*% t(dm)[grep("3.0-3.4", colnames(rnd)), ]
tip_35  <- rnd[, grep("3.5-3.9", colnames(rnd))]   %*% t(dm)[grep("3.5-3.9", colnames(rnd)), ]
tip_4   <- rnd[, grep("4.0-4.4", colnames(rnd))]   %*% t(dm)[grep("4.0-4.4", colnames(rnd)), ]
tip_45  <- matrix(0, ncol = n, nrow = N)  #Reference Interval

#Generate Precision Predictor
precision_lp <- rnd[, grep("gamma", colnames(rnd))] %*% t(dp)

#Convert to means
mu_output <- array(NA, dim = c(n, 4, N),
                   dimnames = list(1:n,
                                   c("4.5-4.9", "4.0-4.4", "3.5-3.9", "3.0-3.4"),
                                   1:N))

for (i in 1:N) {
  denom <- exp(tip_45[i, ]) + exp(tip_4[i, ]) + exp(tip_35[i, ]) + exp(tip_3[i, ])
  
  mu_output[, "4.0-4.4", i] <- exp(tip_4[i, ]) / denom
  mu_output[, "3.5-3.9", i] <- exp(tip_35[i, ]) / denom
  mu_output[, "3.0-3.4", i] <- exp(tip_3[i, ]) / denom
  mu_output[, "4.5-4.9", i] <- 1 - rowSums(mu_output[, c("4.0-4.4", "3.5-3.9", "3.0-3.4"), i])
}

#Take draws from Dirichlet Distribution
alphas <- lapply(dimnames(mu_output)[[2]], function(name) {
  mu_output[, name, ] * t(exp(precision_lp))
})
names(alphas) <- dimnames(mu_output)[[2]]

#Array to store Draws
output <- array(NA, dim = c(n, 4, N), dimnames = dimnames(mu_output))

for (i in 1:N) {
  test <- sapply(alphas, function(a) a[, i])
  output[, , i] <- t(apply(test, 1, function(a) MCMCpack::rdirichlet(1, a)))
}

#Calculate means and confidence intervals across draws
mean.quant <- function(x, probs = c(0.055, 0.945)) {
  c(mean = mean(x), quantile(x, probs = probs))
}

#Grab across all
quant_list <- lapply(1:4, function(j) {
  apply(output[, j, ], 1, mean.quant)
})

names(quant_list) <- dimnames(mu_output)[[2]]

#Compile all ze data
df_mod_tips <- do.call(rbind, lapply(names(quant_list), function(name) {
  q <- t(quant_list[[name]])
  df <- data.frame(df_sp_tis_covar[, c("year", "Region")],
                   trophic_interval = name,
                   mean = q[, "mean"],
                   lower = q[, "5.5%"],
                   upper = q[, "94.5%"])
  return(df)
}))

df_mod_tips_mean <- df_mod_tips %>%
  group_by(Region, trophic_interval) %>%
  summarise(mean_tips = mean(mean),
            sd_tips = sd(mean),
            se_tips = standard_error(mean))

##### Figure 4 #####
df_mod_tips_mean$Region <- factor(df_mod_tips_mean$Region, levels=c('NFLS', 'SS', 'NEUS'))

gg_tips_mod_means <- ggplot(df_mod_tips_mean, aes(x = trophic_interval, y = mean_tips, fill = Region)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", size = 0.25) +
  geom_errorbar(aes(ymin = mean_tips - se_tips, ymax = mean_tips + se_tips), 
                position = position_dodge(width = 0.9), width = 0.3, size = 1) +
  scale_fill_viridis_d() +
  labs(x = "Trophic Interval",
       y = "Mean Proportion") +
  scale_y_continuous(breaks = scales::pretty_breaks(), expand = expansion(mult = c(0.0, 0.20))) +
  theme_bw(base_size = 16) +
  theme(legend.position = "none")

gg_tips_mod_means


# ggsave("../Scale_Dependence_Structural_Stability_LMEs/Figures/Figure SX - Model Estimates - Mean Diffs.jpeg", plot = gg_tips_mod_means, width = 10, height = 4, dpi = 300)

##### Precision Model Mean Differences - Top Model #####





#Estimate Precision
precision_lp <- rnd[, grep("gamma", colnames(rnd))] %*% t(dp)
phi_sim <- exp(precision_lp)

#Compute Confidence Intervals
phi_ci <- apply(phi_sim, 2, function(x) {
  c(mean = mean(x), lower = quantile(x, 0.055), upper = quantile(x, 0.955))
})

phi_df <- as.data.frame(t(phi_ci))

df_precision_reg <- data.frame(df_sp_tis_covar[, c("year", "Region")],
                 mean = phi_df[, "mean"],
                 lower = phi_df[, "lower.5.5%"],
                 upper = phi_df[, "upper.95.5%"])

df_precision_reg_means <- df_precision_reg %>%
  group_by(Region) %>%
  summarise(mean_precis = mean(mean),
            sd_precis = sd(mean),
            se_precis = standard_error(mean))
  
df_precision_reg_means$Region <- factor(df_precision_reg_means$Region, levels=c('NFLS', 'SS', 'NEUS'))

gg_precis_mod_means <- ggplot(df_precision_reg_means, aes(x = Region, y = mean_precis, fill = Region)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", size = 0.25) +
  geom_errorbar(aes(ymin = mean_precis - se_precis, ymax = mean_precis + se_precis), 
                position = position_dodge(width = 0.9), width = 0.3, size = 1) +
  scale_fill_viridis_d() +
  labs(x = "Region",
       y = "Mean Precision (Temporal Variability)") +
  theme_bw(base_size = 16) +
  scale_y_continuous(breaks = scales::pretty_breaks(), expand = expansion(mult = c(0.0, 0.20))) +
  theme(legend.position = "none")

gg_precis_mod_means

# ggsave("../Scale_Dependence_Structural_Stability_LMEs/Figures/Figure 4B - Model Estimates - Precision Diffs.jpeg", plot = gg_precis_mod_means, width = 5, height = 4, dpi = 300)

mod_final_85conf_int <- confint.DirichletRegModel(mod_final, level = 0.85)
mod_final_89conf_int<- confint.DirichletRegModel(mod_final, level = 0.89)

df_precision_coefs <- data.frame(parameter = c("Commercial Harvest Trend", "NAO"),
                                 estimate = c(mod_final_85conf_int$coefficients$gamma$gamma[2], mod_final_85conf_int$coefficients$gamma$gamma[3]),
                                 lower_85 = c(mod_final_85conf_int$ci[[1]]$gamma[2,1], mod_final_85conf_int$ci[[1]]$gamma[3,1]),
                                 upper_85 = c(mod_final_85conf_int$ci[[1]]$gamma[2,2], mod_final_85conf_int$ci[[1]]$gamma[3,2]),
                                 lower_89 = c(mod_final_89conf_int$ci[[1]]$gamma[2,1], mod_final_89conf_int$ci[[1]]$gamma[3,1]),
                                 upper_89 = c(mod_final_89conf_int$ci[[1]]$gamma[2,2], mod_final_89conf_int$ci[[1]]$gamma[3,2])
)


df_precision_coefs$parameter <- factor(df_precision_coefs$parameter, levels=c('NAO', 'Commercial Harvest Trend'))

gg_precis_coefs <- ggplot(df_precision_coefs, aes(y = parameter, x = estimate)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.4) +
  geom_errorbarh(aes(xmin = lower_89, xmax = upper_89), height = 0, linewidth = 1, color = "black") +
  geom_errorbarh(aes(xmin = lower_85, xmax = upper_85), height = 0, linewidth = 4, color = 'black') +
  geom_point(size = 5, shape = 21, stroke = 0.5, color = "black", fill = "aliceblue") +
  xlim(-0.5, 0.5) +
  theme_bw(base_size = 16) +
  theme(axis.text.y = element_blank())

gg_precis_coefs

# ggsave("../Scale_Dependence_Structural_Stability_LMEs/Figures/Figure 4C - Model Estimates - Precision Coefs.jpeg", plot = gg_precis_coefs, width = 5, height = 4, dpi = 300)


gg_precis_grid <- plot_grid(gg_precis_mod_means, gg_precis_coefs, nrow = 1, align = "hv")
gg_precis_grid

# ggsave("../Scale_Dependence_Structural_Stability_LMEs/Figures/Figure 4B,C - Model Estimates - Precision Diffs, Coefs.jpeg", plot = gg_precis_grid, width = 10, height = 4, dpi = 300)

# gg_fig4_grid <- plot_grid(gg_tips_mod_means, gg_precis_grid, nrow = 2, align = "hv")
# gg_fig4_grid

# ggsave("../Scale_Dependence_Structural_Stability_LMEs/Figures/Figure 4 Panel.jpeg", plot = gg_fig4_grid, width = 10, height = 8, dpi = 300)



##### Precision Over Time - Average Across Top 4 Models #####
#Top 4 Models -  2134, 330, 166, 1478

##### Mod 1 - 2134 #####
means <- unlist(coef(models_interaction[[2134]]))
vc <- vcov(models_interaction[[2134]])
N <- 10000
rnd <- rmvnorm(N, mean = means, sigma = vc)
#Design matrix (precision model)
dp <- models_interaction[[2134]]$Z    

#Generate Precision Predictor
precision_lp_1 <- rnd[, grep("gamma", colnames(rnd))] %*% t(dp)

##### Mod 2 - 330 #####
means <- unlist(coef(models_interaction[[330]]))
vc <- vcov(models_interaction[[330]])
N <- 10000
rnd <- rmvnorm(N, mean = means, sigma = vc)
dp <- models_interaction[[330]]$Z

#Generate Precision Predictor
precision_lp_2 <- rnd[, grep("gamma", colnames(rnd))] %*% t(dp)

##### Mod 3 - 166 #####
means <- unlist(coef(models_interaction[[166]]))
vc <- vcov(models_interaction[[166]])
N <- 10000
rnd <- rmvnorm(N, mean = means, sigma = vc)
dp <- models_interaction[[166]]$Z

#Generate Precision Predictor
precision_lp_3 <- rnd[, grep("gamma", colnames(rnd))] %*% t(dp)

##### Mod 3 - 1478 #####
means <- unlist(coef(models_interaction[[1478]]))
vc <- vcov(models_interaction[[1478]])
N <- 10000
rnd <- rmvnorm(N, mean = means, sigma = vc)
dp <- models_interaction[[1478]]$Z

#Generate Precision Predictor
precision_lp_4 <- rnd[, grep("gamma", colnames(rnd))] %*% t(dp)

#Collect Precision Estimates
phi_1 <- exp(precision_lp_1)
phi_2 <- exp(precision_lp_2)
phi_3 <- exp(precision_lp_3)
phi_4 <- exp(precision_lp_4)

weights <- c(0.18, 0.16, 0.09, 0.05)
weights <- weights / sum(weights)

n_keep <- 5000
phi_all <- list(phi_1, phi_2, phi_3, phi_4)

phi_weighted <- do.call(rbind, lapply(seq_along(phi_all), function(i) {
  n_i <- round(weights[i] * n_keep)
  phi_all[[i]][sample(1:nrow(phi_all[[i]]), n_i, replace = TRUE), ]
}))

phi_ci <- apply(phi_weighted, 2, function(x) {
  c(mean = mean(x), lower = quantile(x, 0.055), upper = quantile(x, 0.955))
})

phi_df <- as.data.frame(t(phi_ci))
df_precision_reg <- cbind(df_sp_tis_covar[, c("year", "Region")], phi_df)

df_precision_reg <- df_precision_reg %>%
  rename(lower = `lower.5.5%`,
         upper = `upper.95.5%`)

##### Figure 5 #####
#Regional Grounfish Biomass Collapses
df_gfb_windows <- tribble(
  ~Region, ~xintercept,
  "NFLS",  1990,
  "NFLS",  1994,
  "SS",    1992,
  "SS",    1998,
  "NEUS",  1981,
  "NEUS",  1991,
)

df_mod_tips$Region <- factor(df_mod_tips$Region, levels=c('NFLS', 'SS', 'NEUS'))
df_gfb_windows$Region <- factor(df_gfb_windows$Region, levels=c('NFLS', 'SS', 'NEUS'))

gg_pred_mean_tips <- ggplot(df_mod_tips, aes(x = year, y = mean, color = trophic_interval)) +
  geom_vline(data = df_gfb_windows, aes(xintercept = xintercept), linetype = "dashed", alpha = 0.8) +
  geom_ribbon(aes(ymin = lower, ymax = upper, color = NA, fill = trophic_interval), alpha = 0.25) +
  geom_line(linewidth = 2, alpha = 0.8) +
  facet_wrap(~ trophic_interval + Region, nrow = 4) +
  theme_bw(base_size = 14) +
  scale_color_manual(values = c("#4245c4", "#f78c2f",  "#649f4d", "#b61790")) +
  scale_fill_manual(values = c("#4245c4", "#f78c2f",  "#649f4d", "#b61790")) +
  scale_x_continuous(limits = c(1970, 2005), breaks = scales::pretty_breaks(), expand = expansion(mult = c(0.10, 0.10))) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  labs(color = "Region",
       y = "Mean Trophic Interval Proportions",
       x = "Year") +
  theme(legend.position = "none")

gg_pred_mean_tips

df_precision_reg$Region <- factor(df_precision_reg$Region, levels=c('NFLS', 'SS', 'NEUS'))
df_gfb_windows$Region <- factor(df_gfb_windows$Region, levels=c('NFLS', 'SS', 'NEUS'))

gg_pred_precis <- ggplot(df_precision_reg, aes(x = year, y = mean, color = Region)) +
  geom_vline(data = df_gfb_windows, aes(xintercept = xintercept), linetype = "dashed", alpha = 0.8) +
  geom_ribbon(aes(ymin = lower, ymax = upper, color = NA), alpha = 0.25) +
  geom_line(linewidth = 2, alpha = 0.8) +
  facet_wrap(~Region, nrow = 1, scale = "free_y") +
  theme_bw(base_size = 16) +
  scale_color_viridis_d() +
  scale_x_continuous(limits = c(1970, 2005), breaks = scales::pretty_breaks()
                     ) +
  scale_y_continuous(breaks = scales::pretty_breaks(), 
                     # expand = expansion(mult = c(0.25, 0.25))
                     ) +
  labs(color = "Region",
       y = expression(phi),
       x = "Year") +
  theme(legend.position = "none")

gg_pred_precis

# ggsave("../Scale_Dependence_Structural_Stability_LMEs/Figures/Figure 5 - Model Estimates - Precision.jpeg", plot = gg_pred_precis, width = 10, height = 4, dpi = 300)

# ggsave("../Scale_Dependence_Structural_Stability_LMEs/Figures/Figure 5 - Model Estimates - Precision - Stack.jpeg", plot = gg_pred_precis, width = 4, height = 9, dpi = 300)

# ggsave("../Scale_Dependence_Structural_Stability_LMEs/Figures/Figure SX - Model Estimates - Mean.jpeg", plot = gg_pred_mean_tips, width = 10, height = 10, dpi = 300)

