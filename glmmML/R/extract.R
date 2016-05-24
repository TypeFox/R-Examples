extractAIC.glmmML <- function(fit, scale, k = 2, ...) 
{
    if (k != 2) warning("Only k = 2 is implemented")
    edf <- length(fit$coefficients) + 1
    c(edf, fit$aic)
}

extractAIC.glmmboot <- function(fit, scale, k = 2, ...) 
{
    if (k != 2) warning("Only k = 2 is implemented")
    edf <- length(fit$coefficients) + length(fit$frail)
    c(edf, fit$aic)
}

nobs.glmmML <- function(object, ...){
    object$n
}

nobs.glmmboot <- function(object, ...){
    object$n
}
