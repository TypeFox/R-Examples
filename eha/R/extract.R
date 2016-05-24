extractAIC.coxreg <- function(fit, scale, k = 2, ...) 
{
    edf <- sum(fit$df)
    loglik <- fit$loglik[length(fit$loglik)]
    c(edf, -2 * loglik + k * edf)
}

extractAIC.phreg <- function(fit, scale, k = 2, ...) 
{
    edf <- sum(fit$df)
    loglik <- fit$loglik[length(fit$loglik)]
    c(edf, -2 * loglik + k * edf)
}

extractAIC.aftreg <- function(fit, scale, k = 2, ...) 
{
    edf <- sum(fit$df)
    loglik <- fit$loglik[length(fit$loglik)]
    c(edf, -2 * loglik + k * edf)
}

nobs.coxreg <- function(object, ...){
    object$n
}

nobs.phreg <- function(object, ...){
    object$n
}
nobs.aftreg <- function(object, ...){
    object$n
}

extractAIC.glmmML <- function(fit, scale, k = 2, ...) {
    if (k != 2) warning("Only k = 2 is implemented")
    edf <- length(fit$coefficients) + 1
    c(edf, fit$aic)
}

extractAIC.glmmboot <- function(fit, scale, k = 2, ...) {
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
