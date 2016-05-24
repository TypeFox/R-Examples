################################
#
# THE CONTINUOUSLY RANKED PROBABILITY SCORE FOR ENSEMBLE FORECASTS
#
# ens ... ensemble values (matrix of dimension N*K)
# obs ... observations (vector of length N)
#
################################
EnsCrps <- function(ens, obs) {

  # pre-processing
  l <- Preprocess(ens=ens, obs=obs)
  ens <- l[["ens"]]
  obs <- l[["obs"]]

  N <- length(obs)
  K <- apply(ens, 1, function(x) length(x[!is.na(x)]))
  K[K==0] <- NA

  crps <- sapply(1:N, function(i) 
            mean(abs(ens[i,] - obs[i]), na.rm=TRUE) - 
            sum(dist(ens[i,]), na.rm=TRUE) / (K[i]^2))
  crps[is.nan(crps)] <- NA
  
  return(crps)
}

