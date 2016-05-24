################################
#
# THE FAIR BRIER SCORE
#
# ens ... ensemble matrix (N*K)
# obs ... observation vector (N)
# tau ... exceedance threshold that defines the event
#
# reference: Ferro (2013) "Fair scores for ensemble forecasts"
#
################################
FairBrier <- function(ens, obs, tau=0.5) {

  # pre-processing
  l <- Preprocess(ens=ens, obs=obs)
  ens <- l[["ens"]]
  obs <- l[["obs"]]

  # calculate brier scores #
  ##########################

  # count number of non-NA members per forecast
  K <- apply(ens, 1, function(x) length(x[!is.na(x)]))

  # K=0 or K=1 gives fair brier score NA
  K[K==0 | K==1] <- NA

  # count number of ensemble members > tau
  i <- rowSums(ens > tau, na.rm=TRUE)

  # calculate "observation > tau" indicator
  j <- 1 * (obs > tau)

  # calculate brier score
  fb <- (j - i / K) ^ 2 - i * (K - i) / K / K / (K - 1)
  return(fb)
}

