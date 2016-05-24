###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2012-2014
#
mvSim <- 
function(nsim=1, mu, Sigma, posdeftol=sqrt(.Machine$double.eps), drop=TRUE) {
#
################################################################################
#
  # DIMENSIONS
  k <- length(mu)
  eigen <- eigen(Sigma,symmetric=TRUE)
  # EIGENVALUE DECOMPOSITION
  ev <- eigen$values
#
  # CHECK POSITIVE DEFINITENESS
  if(any(ev < -posdeftol*abs(max(ev)))) stop("lack of positive definiteness")
#
  # SIMULATE
  X <- matrix(rnorm(k*nsim),nsim)
  sim <- drop(mu) + eigen$vectors %*% diag(sqrt(pmax(ev,0)),k) %*% t(X)
#
  if(drop) drop(t(sim)) else t(sim)
}
