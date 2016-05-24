print.mexhaz <- function(x, ...){

    se <- sqrt(diag(x$vcov))
    tval <- coef(x)/se
    df <- x$n.obs-x$n.par
    n.miss <- x$n.obs.tot-x$n.obs

    TAB <- cbind(Estimate = coef(x),
                 StdErr = se,
                 t.value = tval,
                 p.value = 2*pt(-abs(tval),df=df))

    cat("Call:\n")
    print(x$call)
    cat("\nCoefficients:\n")
    printCoefmat(TAB,P.values=TRUE,has.Pvalue=TRUE)
    if(n.miss>0){
        Miss <- paste(n.miss," observations were deleted due to missingness\n",sep="")
    }
    else {
        Miss <- NULL
    }
    if (x$n.time.0>0){
        Time0 <- paste(x$n.time.0," observations had a follow-up time equal to 0 (replaced by 1/730.5)\n",sep="")
    }
    else {
        Time0 <- NULL
    }
    cat(paste("\nlog-likelihood: ",round(x$loglik,4)," (for ",x$n.par," degree(s) of freedom)\n",sep=""))
    cat(paste("\nnumber of observations: ",x$n.obs,", number of events: ",x$n.events,"\n",sep=""))
    cat(Miss)
    cat(Time0)
}
