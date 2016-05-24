logLik.metaplus <- function(object, ...)
{
    logLik(object$fittedmodel)
}

BIC.metaplus <-
function (object, ...) 
{
    if (!is.element("metaplus", class(object))) 
        stop("Argument 'object' must be an object of class \"metaplus\".")
    BIC(object$fittedmodel)
}

AIC.metaplus <-
function(object,...) {
    if (!inherits(object, "metaplus"))
        stop("Use only with 'metaplus' objects.\n")
    AIC(object$fittedmodel)
}
