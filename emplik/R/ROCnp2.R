#### Oct. 2011.
#### For two samples of right censored data.
#### Both samples use nonparametric likelihood = EL. 
#### We test the Ho: ROC curve  R(t0) = b0
#### Or testing the (1-b0)th quantile of sample one equals to
#### the (1-t0)th quantile of sample two.
#### The original version ROCnp try to minimize using the R
#### function optimize( ). Due the the fact that the functions 
#### are piecewise constants, the optimize() may find a false min.
#### This version does an exhaustive search. Slow but find the min.
#### input: t1; right censored times, sample 1.
####        d1; censoring status, d1=1 means uncensored.

ROCnp2 <- function(t1, d1, t2, d2, b0, t0) {
    if (length(b0) != 1) 
        stop("check length of b0")
    if (length(t0) != 1) 
        stop("check length of t0")
    if (b0 >= 1 | b0 <= 0) 
        stop("check the value of b0")
    if (t0 >= 1 | t0 <= 0) 
        stop("check the value of t0")
    tempnp2 <- WKM(x = t2, d = d2)
    place2 <- sum(tempnp2$surv >= t0)
    c2 <- tempnp2$times[place2]
    tempnp1 <- WKM(x = t1, d = d1)
    place1 <- sum(tempnp1$surv > b0)
    c1 <- tempnp1$times[place1]
    if (c2 <= c1) 
        c1 <- tempnp1$times[place1 + 1]
    llr <- function(const, t1, d1, t2, d2, b0, t0) {
        npllik1 <- el.cen.EM2(x = t1, d = d1, fun = function(x, 
            theta) {
            as.numeric(x <= theta)
        }, mu = 1 - b0, theta = const)$"-2LLR"
        npllik2 <- el.cen.EM2(x = t2, d = d2, fun = function(x, 
            theta) {
            as.numeric(x <= theta)
        }, mu = 1 - t0, theta = const)$"-2LLR"
        return(npllik1 + npllik2)
    }
    lower <- min(c1, c2)
    upper <- max(c1, c2)
    timesALL <- c(tempnp2$times, tempnp1$times)
    midpts <- timesALL[(timesALL>lower) & (timesALL<upper)]
    midpts <- sort(midpts)
    midpts <- (midpts[-1] + midpts[-length(midpts)])/2
    k <- length(midpts)
    val2 <- rep(-9, k)
    for(i in 1:k) val2[i] <- llr(const=midpts[i], t1=t1, d1=d1, t2=t2, d2=d2, b0=b0, t0=t0)
    val <- min(val2)
    cstar <- midpts[ which(val2==val) ] 
    list(`-2LLR` = val, cstar = cstar)
}
