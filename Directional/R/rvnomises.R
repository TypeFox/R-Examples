
################################
#### Simulating from a von Mises distribution using the algorithm for
#### the von Mises-Fisher distribution
#### Tsagris Michail 2/2015
#### mtsagris@yahoo.gr
#### References: Wood ATA (1994)
#### Simulation of the von Mises Fisher distribution (Communications in Statistics-Simulation) and
#### Inderjit S. Dhillon and Suvrit Sra (2003)
#### Modeling Data using Directional Distributions (Technical report, The University of Texas at Austin)
################################

rvonmises <- function(n, m, k, rads = TRUE) {
  ## n is the sample size
  ## m is the mean angle expressed in radians or degrees
  ## k is the concentration parameter
  if (rads == F) m <- m/180 * pi  ## turn the degrees into radians
  mu <- c(cos(m), sin(m))
  if (k > 0) {  ## draw from a von Mises distribution
    x <- rvmf(n, mu, k)  ## sample from the von Mises distribution
    ## x is expressed in cartesian coordinates. We have to transform
    ## the values into degrees or radians
    u <- (atan(x[, 2]/x[, 1]) + pi * I(x[, 1] < 0)) %% (2 * pi)  ## u is in radians
  } else u <- runif(n, 0, 2 * pi)  ## draw from a von Mises distribution
  if (rads == F) u <- u * pi/180  ## should the data be in degrees?
  u
}

