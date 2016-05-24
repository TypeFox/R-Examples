################################
#
# THE CONTINUOUSLY RANKED PROBABILITY SCORE FOR GAUSSIAN FORECASTS
#
# mean ... means of the forecast distribution, vector of length N
# sd ... standard deviations of the forecast distributions, vector of lenght N
# obs ... observations, vector of length N
#
# references: * Gneiting, Raftery (2007) "Probabilistic forecasts,
#               calibration and sharpness"
#
################################
GaussCrps <- function(mean, sd, obs) {

  # transform data frames and matrices to vectors
  mean <- c(as.matrix(mean))
  sd <- c(as.matrix(sd))
  obs <- c(as.matrix(obs))

  #
  N <- length(obs)

  # if mean or sd are of length one, expand them to length N
  if (length(mean) == 1) {
    mean <- rep(mean, N)
  }
  if (length(sd) == 1) {
    sd <- rep(sd, N)
  }

  # check if all vectors are of equal length
  len <- c(length(mean), length(sd), length(obs))
  if (length(unique(len)) > 1) {
    stop("vectors have different lengths")
  }

  # initialize crps vector, will be NA whenever sd < 0
  crps <- rep(NA, N)

  # for sd = 0, the crps reduces to the absolute difference
  i.zero <- which(sd == 0)
  crps[i.zero] <- abs(mean[i.zero] - obs[i.zero])

  # crps for sd > 0
  i.pos <- which(sd > 0)
  z <- (obs - mean) / sd
  crps[i.pos] <- (sd * (z * (2 * pnorm(z) - 1) + 
                  2 * dnorm(z) - 1/sqrt(pi)))[i.pos]

  # set all bad values to NA
  # bad values are NA, Inf, NaN, sd<0
  bad <- (is.na(crps) | is.infinite(crps) | is.nan(crps) | sd < 0)
  crps[bad] <- NA

  # return crps vector
  return(crps)
}


