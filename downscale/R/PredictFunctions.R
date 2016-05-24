################################################################################
# 
# PredictFunctions.R
# Version 1.2
# 28/01/2015
#
# Updates:
#   30/01/2015: Thomas model added
#   28/01/2015: Altered gamma function in Finite Negative Binomial model
#
# Functions for downscaled prediction of area of occupancy for a given 
# grain size (A) given the parameters for that model. Model parameters may
# be the outputs from model optimisation to coarse-scale data
#
# This file contains the 8 functions that are simple geometric extrapolations of
# the area-occupancy relationship at coarse grain sizes:
#   Nachman   Nachman model
#   PL        Power Law model
#   Logis     Logistic model
#   Poisson   Poisson model
#   NB        Negative binomial model
#   GNB       Generalised negative binomial model
#   INB       Improved negative binomial model
#   FNB       Finite negative binomial model
#   Thomas    Thomas model
#
# requires(Rmpfr)
#   # Note: The Finite Negative Binomial model requires multiple precision 
#     floating-point numbers in package Rmpfr for calculation (see description)
#
################################################################################

### Nachman
PredictNachman <- function(par, area) {
  # Predicts area of occupancy for grain size A using the Nachman model
  #
  # Args:
  #   par: dataframe containing parameters C and z of the Nachman model
  #   area: Grain size (km2) to be predicted
  AOO <- log(1 - exp(-par$C * area ^ par$z))
  return(AOO)
}

### Power Law
PredictPL <- function(par, area) {
  # Predicts area of occupancy for grain size A using the Power law model
  #
  # Args:
  #   par: dataframe containing parameters C and z of the Power law model
  #   area: Grain size (km2) to be predicted
  AOO <- log(par$C * area ^ par$z)
  return(AOO)
}

### Logistic
PredictLogis <- function(par, area) {
  # Predicts area of occupancy for grain size A using the Logistic model
  #
  # Args:
  #   par: dataframe containing parameters C and z of the Logistic model
  #   area: Grain size (km2) to be predicted
  AOO <- log((par$C * (area ^ par$z)) / (1 + (par$C * (area ^ par$z))))
  return(AOO)
}

### Poisson
PredictPoisson <- function(par, area) {
  # Predicts area of occupancy for grain size A using the Poisson model
  #
  # Args:
  #   par: dataframe containing parameter lambda of the Poisson model
  #   area: Grain size (km2) to be predicted
  AOO <- log(1 - (exp(-par$lambda * area)))
  return(AOO)
}

### Negative binomial model
PredictNB <- function(par, area) {
  # Predicts area of occupancy for grain size A using the Negative 
  # Binomial model
  #
  # Args:
  #   par: dataframe containing parameters C and k of the Negative
  #        Binomial model
  #   area: Grain size (km2) to be predicted
  AOO <- log(1 - (1 + (par$C * area) / par$k) ^ -par$k)
  return(AOO)
}

### Generalised negative binomial model
PredictGNB <- function(par, area) {
  # Predicts area of occupancy for grain size A using the Generalised Negative 
  # Binomial model
  #
  # Args:
  #   par: dataframe containing parameters C, z and k of the Generalised 
  #        Negative Binomial model
  #   area: Grain size (km2) to be predicted
  AOO <- log(1 - (1 + (par$C * area ^ par$z) / par$k) ^ -par$k)
  return(AOO)
}

### Improved negative binomial model
PredictINB <- function(par, area) {
  # Predicts area of occupancy for grain size A using the Improved Negative 
  # Binomial model
  #
  # Args:
  #   par: dataframe containing parameters C and b of the Improved 
  #        Negative Binomial model
  #   area: Grain size (km2) to be predicted
  AOO <- log(1 - ((par$C * area ^ (par$b - 1))^
                ((par$r * area) / (1 - par$C * area ^ (par$b - 1)))))
  return(AOO)
}

### Finite negative binomial model
PredictFNB <- function(par, area, extent){
  # Predicts area of occupancy for grain size A using the Finite Negative 
  # Binomial model. The function multiplies many gamma functions and so
  # integers may become larger than possible in R. Therefore  we use multiple
  # precision floatinf point numbers (the 'mpfr' function in package 'Rmpfr')
  # is used to make calculations possible.
  #
  # Args:
  #   par: dataframe containing parameters W and k of the Finite
  #        Negative Binomial model
  #   area: Grain size (km2) to be predicted
  #   extent: Total area (km2)
  gamma1 <- par$W + ((extent * par$k) / area) - par$k
  gamma2 <- (extent * par$k) / area
  gamma3 <- par$W + ((extent * par$k) / area)
  gamma4 <- ((extent * par$k) / area) - par$k
  AOO <- suppressWarnings(as.numeric(log(1 - (
    (gamma(Rmpfr::mpfr(gamma1, 64)) * gamma(Rmpfr::mpfr(gamma2, 64))) /
      (gamma(Rmpfr::mpfr(gamma3, 64)) * gamma(Rmpfr::mpfr(gamma4, 64)))))))
  return(AOO)
}

### Thomas model
PredictThomas <- function(par, area, extent, tolerance = 1e-6){
  # Predicts area of occupancy for grain size A using the Thomas model.
  #
  # Args:
  #   par: dataframe containing parameters rho, mu and sigma of the Thomas model
  #   area: Grain size (km2) to be predicted
  #   extent: Total area (km2)
  #   tolerance: tolerance of the integration. The smaller the number the
  #     greater the accuracy but longer the processing time
  AOO <- sapply(1:length(area), 
                function(i) log(1 - exp(-par$rho * espA(par = par,
                                                        area = area[i],
                                                        extent = extent,
                                                        tolerance = tolerance)
                                        )))
  return(AOO)
}
