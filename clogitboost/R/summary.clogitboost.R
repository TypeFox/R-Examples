
summary.clogitboost <- function(object, ...){
  res <- list(call = object$call, infscore = object$infscore, loglike=object$loglike)
  res
}
