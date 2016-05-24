summary.mexhaz <- function(object, ...){
    coef <- coef(object)
    se <- sqrt(diag(object$vcov))
    tval <- coef(object)/se
    df <- object$n.obs-object$n.par
    n.miss <- object$n.obs.tot-object$n.obs

    TAB <- cbind(Estimate = coef,
                 StdErr = se,
                 t.value = tval,
                 p.value = 2*pt(-abs(tval),df=df))

    HRTab <- NULL
    ishr <- 0
    IdxPH <- which(names(coef)%in%object$names.ph)
    if (length(IdxPH)>0){
        PH <- coef[IdxPH]
        sePH <- se[IdxPH]
        HR <- exp(PH)
        ci.lower <- exp(PH+qt(0.025,df=df)*sePH)
        ci.upper <- exp(PH+qt(0.975,df=df)*sePH)
        HRTab <- cbind(Coef=round(PH,4),HR=round(HR,4),CI.lower=round(ci.lower,4),CI.upper=round(ci.upper,4))
        ishr <- 1
    }

    cat("Call:\n")
    print(object$call)
    cat("\nCoefficients:\n")
    printCoefmat(TAB,P.values=TRUE,has.Pvalue=TRUE)
    if (ishr){
        cat("\nHazard ratios (for proportional effect variables):\n")
        print(HRTab)
    }
    if(n.miss>0){
        Miss <- paste(n.miss," observations were deleted due to missingness\n",sep="")
    }
    else {
        Miss <- NULL
    }
    if (object$n.time.0>0){
        Time0 <- paste(object$n.time.0," observations had a follow-up time equal to 0 (replaced by 1/730.5)\n",sep="")
    }
    else {
        Time0 <- NULL
    }

    cat(paste("\nlog-likelihood: ",round(object$loglik,4)," (for ",object$n.par," degree(s) of freedom)\n",sep=""))
    cat(paste("\nnumber of observations: ",object$n.obs,", number of events: ",object$n.events,"\n",sep=""))
    cat(Miss)
    cat(Time0)
}
