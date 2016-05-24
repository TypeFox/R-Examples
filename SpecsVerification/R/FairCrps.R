################################
#
# THE FAIR CONTINUOUSLY RANKED PROBABILITY SCORE
#
# ens ... ensemble values (matrix of dimension N*K)
# obs ... observations (vector of length N)
#
# references: * Gneiting, Raftery (2007) "Probabilistic forecasts,
#               calibration and sharpness"
#             * Ferro, Richardson, Weigel (2008) "On the effect 
#               of ensemble size on the discrete and continuous
#               ranked probability scores"
#
################################
FairCrps <- function(ens, obs) {

  # pre-processing
  l <- Preprocess(ens=ens, obs=obs)
  ens <- l[["ens"]]
  obs <- l[["obs"]]

  N <- length(obs)
  K <- apply(ens, 1, function(x) length(x[!is.na(x)]))
  K[K==0 | K==1] <- NA

  crps <- sapply(1:N, function(i) 
            mean(abs(ens[i,] - obs[i]), na.rm=TRUE) - 
            sum(dist(ens[i,]), na.rm=TRUE) / 
            (K[i] * (K[i]-1)))
  crps[is.nan(crps)] <- NA
  
  return(crps)
}

