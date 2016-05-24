##' @export
logLik.biprobit <- function(object,indiv=FALSE,...) {
  if (indiv) return(object$logLik)
  n <- sum(object$N[1])
  p <- length(coef(object))
  loglik <- sum(object$logLik)
  attr(loglik, "nall") <- n
  attr(loglik, "nobs") <- n
  attr(loglik, "df") <- p
  class(loglik) <- "logLik"        
  return(loglik)
}
