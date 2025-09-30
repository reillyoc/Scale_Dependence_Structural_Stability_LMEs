# Functions to Calculate Different Variability Metrics

# Author(s): Reilly O'Connor
# Version: 2024-06-18

# Pkgs
library(tidyverse)
library(zoo)

##### Code #####

#function to calculate standard error
standard_error <- function(x) {
  sd(x) / sqrt(length(x))
}

#Coefficient of Variation (CV)
cv <- function(x) {
  if (length(x) < 2) return(NA)
  CV <- sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
  return(CV)
}

k_cv <- function(x) {
  if (length(x) < 2) return(NA)
  CV <- sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
  k_CV <- sqrt( (CV^2) /(1 + (CV^2)))
  return(k_CV)
}

#Proportional Variability (PV)
p_cv <- function(x) {
  n <- length(x)
  if (n < 2) return(NA)
  z_values <- outer(x, x, function(a, b) 1 - pmin(a, b) / pmax(a, b))
  P_CV <- 2 * sum(z_values[upper.tri(z_values)], na.rm = TRUE) / (n * (n - 1))
  return(P_CV)
}

# Consecutive Disparity Index (D)
disp_cv <- function(x) {
  n <- length(x)
  if (n < 2) return(NA)
  D_CV <- sum(abs(log(x[-1] / x[-n])), na.rm = TRUE) / (n - 1)
  return(D_CV)
}

acf_one <- function(x) {
  if (length(unique(x)) > 1) {
    return(acf(x, lag.max = 1, plot = FALSE)$acf[2])
  } else {
    return(NA_real_)  # Avoid error if data is constant
  }
}

CalcZeroInfGeomDens = function(x){
  prob_obs = sum(x > 0)/length(x)
  geom_mean_dens = ifelse(prob_obs > 0, exp(mean(log(x[x > 0]))), 0)
  return(geom_mean_dens * prob_obs)
}



pyper_peterson_corr <- function(x, y, max_lag = 1, alpha = 0.10, ac_est = TRUE) {
  stopifnot(length(x) == length(y))
  n <- length(x)
  
  # Standardize x and y
  x <- scale(x, center = TRUE, scale = TRUE)[, 1]
  y <- scale(y, center = TRUE, scale = TRUE)[, 1]
  
  # Compute autocorrelations (exclude lag 0)
  acf_x <- acf(x, lag.max = max_lag, plot = FALSE)$acf[-1]
  acf_y <- acf(y, lag.max = max_lag, plot = FALSE)$acf[-1]
  
  # === Apply Equation 7 correction factor ===
  if (ac_est) {
    correction_factors <- sapply(1:max_lag, function(j) n / (n - j))
    acf_x <- acf_x * correction_factors
    acf_y <- acf_y * correction_factors
  }
  
  # Compute weighted sum for neff denominator
  weighted_sum <- sum(sapply(1:max_lag, function(j) {
    ((n - j) / n) * acf_x[j] * acf_y[j]
  }))
  
  # Effective sample size (N*)
  neff_inv <- (1 / n) + (2 / n) * weighted_sum
  neff <- 1 / neff_inv
  neff <- max(neff, 4)  # Prevent overcorrection
  
  # Pearson correlation and confidence intervals
  r <- cor(x, y, method = "pearson")
  z <- 0.5 * log((1 + r) / (1 - r))  # Fisher z-transform
  se <- 1 / sqrt(neff - 3)
  tcrit <- qt(1 - alpha / 2, df = neff - 2)
  lower <- tanh(z - tcrit * se)
  upper <- tanh(z + tcrit * se)
  
  tibble(
    r = r,
    lower = lower,
    upper = upper,
    neff = neff
  )
}

