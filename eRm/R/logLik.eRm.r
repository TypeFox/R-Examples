logLik.eRm <- function(object,...)
{
#object of class eRm
  # val <- object$loglik
  # attr(val, "df") <- object$npar
  val <- list(loglik = object$loglik, df =  object$npar) # rh 26-03-2010
  class(val) <- "logLik.eRm"
  val
}


