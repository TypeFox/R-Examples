##extract log-likelihood of model
##generic
extractLL <- function(mod, ...) {
  UseMethod("extractLL", mod)
}


##generic
extractLL.default <- function(mod, ...) {
stop("\nFunction not yet defined for this object class\n")
}

  
##methods
##coxme objects
extractLL.coxme <- function(mod, type = "Integrated", ...) {
  ##fixed effects
  fixed.K <- length(fixef(mod))
  ##random effects
  random.K <- length(ranef(mod))
  df <- fixed.K + random.K
  LL <- mod$loglik[type]
  attr(LL, "df") <- df
  return(LL)
}


##coxph objects
extractLL.coxph <- function(mod, ...) {
  coefs <- coef(mod)
  if(is.null(coefs)) {
      ncoefs <- 0
      LL <- mod$loglik[1] #when null model, only 1 log-likelihood value
    } else {
      ncoefs <- length(coefs)
      LL <- mod$loglik[2] #second value is the logLik at the solution
    }
  attr(LL, "df") <- ncoefs
  return(LL)
}


##lmekin objects
extractLL.lmekin <- function(mod, ...) {
  LL <- mod$loglik
  #K = fixed + random + residual variance
  fixed.K <- length(fixef(mod))
  random.K <- length(ranef(mod))
  df <- fixed.K + random.K + 1
  attr(LL, "df") <- df
  return(LL)
}


##maxlikeFit objects
extractLL.maxlikeFit <- function(mod, ...) {
  LL <- logLik(mod)
  df <- length(coef(mod))
  attr(LL, "df") <- df
  return(LL)
}


##unmarkedFit objects
extractLL.unmarkedFit <- function(mod, ...) {
  LL <- -1*mod@negLogLike
  df <- length(mod@opt$par)
  attr(LL, "df") <- df
  return(LL)
}


##vglm objects
extractLL.vglm <- function(mod, ...) {
  LL <- logLik(mod)
  df <- length(coef(mod))
  if(identical(mod@family@vfamily, "gaussianff")) {df <- df + 1}
  attr(LL, "df") <- df
  return(LL)
}
