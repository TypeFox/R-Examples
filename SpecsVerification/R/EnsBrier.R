################################
#
# THE BRIER SCORE FOR ENSEMBLE FORECASTS
#
# ens ... ensemble matrix (N*K)
# obs ... observation vector (N)
# tau ... exceedance threshold that defines the event
#
################################
EnsBrier <- function(ens, obs, tau=0.5) {

  # pre-processing
  l <- Preprocess(ens=ens, obs=obs)
  ens <- l[["ens"]]
  obs <- l[["obs"]]

  # calculate brier scores #
  ##########################
  # count number of non-NA members per forecast
  K <- apply(ens, 1, function(x) length(x[!is.na(x)]))
  K[K==0] <- NA
  # count number of ensemble members > tau
  i <- rowSums(ens > tau, na.rm=TRUE)
  # calculate "observation > tau" indicator
  j <- 1 * (obs > tau)
  # calculate brier score
  br <- (j - i / K) ^ 2 
  return(br)
}

