###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2013-2014
#
mvmetaSim <- 
function(y, S, Psi, sd, cor, nsim=1, seed=NULL, posdeftol) {
#
################################################################################
#
  if(missing(posdeftol)) posdeftol <- sqrt(.Machine$double.eps)
#
  # DEFINE THE SEED (FROM simulate.lm)
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) runif(1)
  if (is.null(seed)) {
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  } else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
#
  # PREPARE AND CHECK y
  if(!is.matrix(y)) y <- as.matrix(y)
  k <- ncol(y)
  if(any(is.na(y))) stop("missing values not allowed in 'y'")
#
  # PREPARE AND CHECK S
  S <- mkS(S,y)
  if(any(is.na(S))) stop("missing values not allowed in 'S'")
#
  # PREPARE AND CHECK Psi, THEN SET THE DIMENSION
  if(missing(Psi)) {
    if(missing(sd)||missing(cor)) stop("'Psi' or 'sd'-'cor' must be provided")
    Psi <- inputcov(sd,cor)
  }
  if(!is.matrix(Psi)) Psi <- xpndMat(Psi)
  Psi <- checkPD(Psi,k,"Psi")
#
  # SAMPLE THE RESPONSES
  # FOR EFFICIENCY, IT SAMPLES SEVERAL OUTCOMES FROM THE SAME MEAN AND
  #   THEN RE-ARRANGE THEM
  sim <- do.call("cbind",lapply(seq(nrow(y)), function(i) {
    mvSim(nsim,y[i,],Sigma=xpndMat(S[i,])+Psi,posdeftol=posdeftol,drop=FALSE)}))
  sim <- lapply(seq(nrow(sim)), function(i) drop(matrix(sim[i,],
    ncol=ncol(y),byrow=T,dimnames=dimnames(y))))
  if(nsim==1) sim <- sim[[1]]
#
  sim
}
