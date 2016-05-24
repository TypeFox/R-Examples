AIC.glmssn <- function(object,...,k=2) {
    if(class(object) != "glmssn") return("Not a glmssn object")
    if(!missing(k)) cat("This argument has no effect\n")
    if(object$args$family == "poisson") return(NA)
    if(object$args$family =="binomial") return(NA)
    cp <- covparms(object)[,3]
    nparmstheta <- length(cp)
    rankX <- object$sampinfo$rankX
    if(object$args$EstMeth == "REML")
        nparms <- nparmstheta
    if(object$args$EstMeth == "ML")
        nparms <- rankX + nparmstheta
    object$estimates$m2LL + k*nparms
}

