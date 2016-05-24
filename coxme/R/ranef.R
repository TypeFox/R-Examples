# Automatically generated from all.nw using noweb
# The objects that do the actual work (not much work)
fixef.coxme <- function(object, ...)
    object$coefficients

fixef.lmekin <- function(object, ...)
    object$coefficients$fixed

ranef.coxme <- function(object, ...)
    object$frail

ranef.lmekin <- function(object, ...)
    object$coefficients$random

VarCorr.coxme <- function(x, ...) 
    x$vcoef
    
VarCorr.lmekin <- function(x, ...) 
    x$vcoef

vcov.coxme <- function(object, ...) {
    nf <- length(fixef(object))
    indx <- seq(length=nf, to=nrow(object$var))
    as.matrix(object$var[indx, indx])
}

vcov.lmekin <- vcov.coxme  
logLik.coxme <- function(object, type=c("penalized", "integrated"), ...) {
    type <- match.arg(type)
    if (type=='penalized') {
        out <- object$loglik[3] + object$penalty
        attr(out, "df") <- object$df[2]
        }
    else {
        out <- object$loglik[2]
        attr(out, "df") <- object$df[1]
        }
    attr(out, "nobs") <- object$n[1]  #number of events
    class(out) <- "logLik"
    out
}

logLik.lmekin <- function(object, ...) {
    out <- object$loglik[2]
    attr(out, "df") <- object$df[1]
    class(out) <- "logLik"
    out
    }
