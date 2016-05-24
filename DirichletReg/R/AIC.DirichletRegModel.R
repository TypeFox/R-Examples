AIC.DirichletRegModel <- function(object, ..., k = 2){
  - 2*object$logLik + k*object$npar
}

nobs.DirichletRegModel <- function(object, ...){
  sum(object$weights)
}

BIC.DirichletRegModel <- function(object, ...){
  - 2*object$logLik + log(nobs(object))*object$npar
}
