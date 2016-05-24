myWKM<-function (x, d, zc = rep(1, length(d)), w = rep(1, length(d))) 
{
    if (any((d != 0) & (d != 1))) 
        stop("d must be 0(right-censored) or 1(uncensored)")
    temp <- myWdataclean3(x, d, zc, w)
    dd <- temp$dd
    ww <- temp$weight
    dd[length(dd)] <- 1
    allrisk <- rev(cumsum(rev(ww)))
    survP <- cumprod(1 - (dd * ww)/allrisk)
    jumps <- -diff(c(1, survP))
    logel <- sum(ww[dd == 1] * log(jumps[dd == 1])) + sum(ww[dd == 
        0] * log(survP[dd == 0]))
    list(times = temp$value, jump = jumps, surv = survP, logel = logel,
        weight=ww)
 }