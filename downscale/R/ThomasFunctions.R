################################################################################
# 
# ThomasFunctions.R
# Version 1.0
# 30/01/2015
#
# Set of functions required for calculations in the Thomas model
#
################################################################################

erf <- function(x) {
  # Function for calculating the error in the isotropic gaussian distribution
  return(2 * pnorm(x * sqrt(2)) - 1)
}

isonormdist <- function(z, area, mu, sigma) {
  # Function for defining the isotropic gaussian distribution
  #
  # Args:
  #   z:  provides the two dimensions (u and v)
  #   area: grain size to be estimated (km2)
  #   mu: parameter of the thomas model
  #   sigma: parameter of the thomas model
  u <- z[1]
  v <- z[2]
  return(exp(-mu * (1 / 4) * 
               ((erf((sqrt(area) - (2 * u)) / (2 * sqrt(2) * sigma))) + 
                  (erf((sqrt(area) + (2 * u)) / (2 * sqrt(2) * sigma)))) * 
               ((erf((sqrt(area) - (2 * v)) / (2 * sqrt(2) * sigma))) +
                  (erf((sqrt(area) + (2 * v)) / (2 * sqrt(2) * sigma))))))  
}

espA <- function(par, area, extent, tolerance = 1e-6) {
  # Function to integrate the 2 dimensional normal distribution
  #
  # Args:
  #   par:  parameters of the Thomas model
  #   area: grain size to be estimated (km2)
  #   extent: total area (km2)
  #   tolerance: tolerance during the integration. A smaller tolerance allows
  #              for greater accuracy but longer processing times
  ext <- extent
  return(ext - cubature::adaptIntegrate(isonormdist,
                                    area = area,
                                    mu = par$mu,
                                    sigma = par$sigma,
                                    tol = tolerance,
                                    lowerLimit = rep(-sqrt(ext) / 2,2),
                                    upperLimit = rep(sqrt(ext) / 2,2)) [[1]])
}

