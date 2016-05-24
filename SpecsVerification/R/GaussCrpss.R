################################
#
# ANALYZE DIFFERENCE IN THE CRPS BETWEEN TWO FORECASTS ISSUED AS NORMAL
# DISTRIBUTIONS FOR THE SAME OBSERVATIONS BY SKILL SCORE 1 - S / S.ref
#
# mean     ... forecast mean (vector of length N)
# sd       ... forecast standard deviation (vector of length N)
# mean.ref ... mean of the reference forecast (vector of length N)
# sd.ref   ... standard deviation of the reference forecast (vector of length N)
# obs      ... observations (vector of length N)
#
################################
GaussCrpss <- function(mean, sd, mean.ref, sd.ref, obs) {

  # transform data frames and matrices to vectors
  mean <- c(as.matrix(mean))
  mean.ref <- c(as.matrix(mean.ref))
  sd <- c(as.matrix(sd))
  sd.ref <- c(as.matrix(sd.ref))
  obs <- c(as.matrix(obs))

  # 
  N <- length(obs)

  # if any forecast vector is of length one, expand them to length N
  if (length(mean) == 1) {
    mean <- rep(mean, N)
  }
  if (length(sd) == 1) {
    sd <- rep(sd, N)
  }
  if (length(mean.ref) == 1) {
    mean.ref <- rep(mean.ref, N)
  }
  if (length(sd.ref) == 1) {
    sd.ref <- rep(sd.ref, N)
  }
  

  # check if all vectors are of equal length
  len <- c(length(mean), length(sd), length(obs), length(mean.ref), length(sd.ref))
  if (length(unique(len)) > 1) {
    stop("vectors have different lengths")
  }

  # calculate crps difference
  crps.ens <- GaussCrps(mean, sd, obs)
  crps.ref <- GaussCrps(mean.ref, sd.ref, obs)

  crpss <- 1 - mean(crps.ens) / mean(crps.ref)
  crpss.sigma <- 1 / sqrt(N) * sqrt( var(crps.ens) / mean(crps.ref)^2 + 
         var(crps.ref) * mean(crps.ens)^2 / mean(crps.ref)^4 - 
         2 * cov(crps.ens, crps.ref) * mean(crps.ens) / mean(crps.ref)^3)

  #return
  list(crpss=crpss, crpss.sigma=crpss.sigma)
}




