logLik.manyglm <- function(object, ...)
{
  logL = (object$two.loglike)/2
  names(logL) = names(object$fitt)
  attributes(logL)$df <- (object$aic[1]+object$two.loglike[1])/2
  attributes(logL)$nobs <- dim(object$fitt)[1]
  class(logL) <- "logLik"
  return( logL )
}

