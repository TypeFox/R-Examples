#----------------------------------------------------------------------------------------
# this a generic function for logLik for gamlss
# it should moved to gamlss extra
logLik.gamlss <- function(object, ...)
{
if (!is.gamlss(object))  stop(paste("This is not an gamlss object", "\n", ""))
  val <- -object$G.deviance/2 #sum(lik)
 attr(val, "nall") <- object$N
    attr(val, "nobs") <- object$noObs
    attr(val, "df") <- object$df.fit
    class(val) <- "logLik"
    val
}             
