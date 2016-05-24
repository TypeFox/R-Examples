se <- function(object, ...) UseMethod("se")

se.tsglm <- function(object, B, parallel=FALSE, ...){
  tsglm.check(object)
  est <- c(coef(object), sigmasq=if(object$distr=="poisson") NULL else object$sigmasq)
  if(missing(B)){
    vcov <- vcov(object)
    var <- diag(vcov)
    stderrors <- c(sqrt(var), sigmasq=if(object$distr=="poisson") NULL else NA)
    result <- list(est=est, se=stderrors, type="normapprox")
  }else{
    stopifnot(B>=2, B%%1==0)
    simfit <- function(seed, fit, ...){
      set.seed(seed)
      ts_sim <- tsglm.sim(fit=fit)$ts
      fit_sim <- tsglm(ts=ts_sim, model=fit$model, xreg=fit$xreg, link=fit$link, distr=fit$distr, score=FALSE, info="none", ...)
      result <- c(coef(fit_sim), sigmasq=if(object$distr=="poisson") NULL else fit_sim$sigmasq)
      return(result)
    }
    seeds <- sample(1e+9, size=B)
    if(parallel){
      Sapply <- function(X, FUN, ...) parSapply(cl=NULL, X=X, FUN=FUN, ...)
    }else{
      Sapply <- sapply
    }
    bootstrap_coefs <- Sapply(seeds, simfit, fit=object, ..., simplify=TRUE)
    if(object$distr!="poisson" && anyNA(bootstrap_coefs["sigmasq",])) warning(paste("The overdispersion coefficient 'sigmasq' could not be estimated\nin", sum(is.na(bootstrap_coefs["sigmasq",])), "of the", B, "replications. It is set to zero for these\nreplications. This might to some extent result in an overestimation\nof its true variability."))
    stderrors <- apply(bootstrap_coefs, 1, sd, na.rm=TRUE)
    result <- list(est=est, se=stderrors, type="bootstrap", B=B)
  }
  return(result)
}
