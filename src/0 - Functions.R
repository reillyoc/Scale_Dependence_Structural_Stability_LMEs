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

#Autocorrelation Lag 1
acf_one <- function(x) {
  if (length(unique(x)) > 1) {
    return(acf(x, lag.max = 1, plot = FALSE)$acf[2])
  } else {
    return(NA_real_)  # Avoid error if data is constant
  }
}

#Function from Pedersen et al 2017, zero-inflated geometric mean density
CalcZeroInfGeomDens = function(x){
  prob_obs = sum(x > 0)/length(x)
  geom_mean_dens = ifelse(prob_obs > 0, exp(mean(log(x[x > 0]))), 0)
  return(geom_mean_dens * prob_obs)
}


