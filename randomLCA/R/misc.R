logLik.randomLCA <- function(object, ...)
{
     val <- object$logLik
    attr(val, "df") <- object$np
    attr(val, "nobs") <- object$nobs
    class(val) <- "logLik"
    val
}

BIC.randomLCA <-
function (object, ...) 
{
    if (!is.element("randomLCA", class(object))) 
        stop("Argument 'object' must be an object of class \"randomLCA\".")
    BIC(logLik(object))
}

AIC.randomLCA <-
function(object,...,k=2) {
    if (!inherits(object, "randomLCA"))
        stop("Use only with 'randomLCA' objects.\n")
    AIC(logLik(object),...,k)
}

AIC3 <-
  function(object) {
    if (!inherits(object, "randomLCA"))
      stop("Use only with 'randomLCA' objects.\n")
    AIC(logLik(object),k=3)
  }
