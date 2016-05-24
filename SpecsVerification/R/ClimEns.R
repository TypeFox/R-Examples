###################################
#
# CREATE A CLIMATOLOGICAL ENSEMBLE FROM A VECTOR OF OBSERVATIONS
#
# INPUT: OBS ... A VECTOR OF LENGTH N
# OUTPUT: ENS.CLIM ... A MATRIX WITH N ROWS AND (N-1) COLUMNS
#
###################################

ClimEns <- function(obs, leave.one.out=TRUE) {

  obs <- Preprocess(obs=obs)[["obs"]]

  if (length(obs) < 2 & leave.one.out == TRUE) {
    stop("Need at least 2 observations to construct leave-one-out ensemble")
  }

  # construct climatological ensemble matrix
  N <- length(obs)
  ens.clim <- t(matrix(rep(obs, N), N, N))
  
  # remove diagonal if desired
  if (leave.one.out) {
    ens.clim <- t(matrix(t(ens.clim)[-seq(1,N^2,N+1)], N-1, N))
  }
  
  return(ens.clim)
}

