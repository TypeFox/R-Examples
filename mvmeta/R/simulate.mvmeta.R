###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2014
#
simulate.mvmeta <- 
function(object, nsim=1, seed=NULL, ...) {
#
################################################################################
#
  # ERROR FOR FUTURE IMPLEMENTATION OF MIKE'S MODEL
  if(!object$method%in%c("fixed","ml","reml","mm","vc"))
    stop("Simulating responses not allowed for estimation method used in model fitting")
#
  # DEFINE THE SEED (FROM simulate.lm)
  if(!exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE)) runif(1)
  if(is.null(seed)) {
    RNGstate <- get(".Random.seed", envir=.GlobalEnv)
  } else {
    R.seed <- get(".Random.seed", envir=.GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind=as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir=.GlobalEnv))
  }
#
  # RECOVER FITTED VALUES AND S (EXCLUDING MISSING VALUES) + Psi
  fit <- as.matrix(na.omit(fitted(object)))
  S <- if(is.null(nas <- attr(fit,"na.action")))
    object$S else object$S[-nas,,drop=FALSE]
  Psi <- if(is.null(object$Psi)) matrix(0,ncol(fit),ncol(fit)) else object$Psi
#
  # FOR EFFICIENCY, IT SAMPLES SEVERAL OUTCOMES FROM THE SAME MEAN AND
  #   THEN RE-ARRANGE THEM
  sim <- do.call("cbind",lapply(seq(nrow(fit)), function(i) {
    mvSim(nsim,fit[i,],Sigma=xpndMat(S[i,])+Psi,drop=FALSE)}))
  sim <- lapply(seq(nrow(sim)), function(i) drop(matrix(sim[i,],
    ncol=ncol(fit),byrow=T,dimnames=dimnames(fit))))
  if(nsim==1) sim <- sim[[1]]
#
  sim
}
