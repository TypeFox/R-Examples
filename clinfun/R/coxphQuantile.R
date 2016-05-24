coxphQuantile <- function(phfit, xrange, p=0.5, whichx=1, otherx=NULL, ...) {
    requireNamespace("survival")
    if (class(phfit) != "coxph") stop("phfit shoud be coxph class object")
    cvtmean <- phfit$means
    loghr <- phfit$coef
    S0 <- survival::survfit(phfit)
    ii <- S0$surv < 1 & S0$surv > 0
    stime <- S0$time[ii]
    ssurv <- S0$surv[ii]
    if (!missing(otherx)) {
        ssurv <- ssurv^(exp(sum(loghr[-whichx]*(otherx-cvtmean[-whichx]))))
    }
    sx <- cvtmean[whichx] + log(log(p)/log(ssurv))/loghr[whichx]
    ii <- which(sx >= xrange[1] & sx <= xrange[2])
    lines(sx[ii], stime[ii], type="S", ...)
    invisible(as.data.frame(list(cvt=sx,time=stime)))
}
