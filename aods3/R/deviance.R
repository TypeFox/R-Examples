deviance.aodml <- function(object, ...) 2 * sum(object$lmax - object$l) # 2 * (logLmax - logL)

deviance.aodql <- function(object, ...) deviance(object$fm, ...)

