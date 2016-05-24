################################
#
# THE FAIR RANKED PROBABILITY SCORE
#
# ens ... ensemble matrix, ens[i,j] = number of ensemble members that predict category j at time i
# obs ... observation matrix, obs[i,j] = 1 if category j is observed at time i, 0 otherwise
#
# reference: Ferro (2013) "Fair scores for ensemble forecasts"
#
################################
FairRps <- function(ens, obs) {

  # preprocess
  l <- Preprocess(ens=ens)
  ens <- l[["ens"]]

  stopifnot(all(dim(ens)==dim(obs)))
  stopifnot(all(rowSums(obs) == 1))

  frps <- 
  sapply(1:nrow(ens), function(i) {
    E <- cumsum(ens[i,])
    O <- cumsum(obs[i,])
    M <- tail(E,1)
    sum((O-E/M)^2 - E * (M - E) / (M*M*(M-1)))
  })
  return(frps)
}

