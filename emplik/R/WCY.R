WCY <- function(x, d, zc = rep(1, length(d)),
         wt = rep(1,length(d)), maxit = 25, error = 1e-09)
{
### Chang and Yang self-consistant est. for doubly censored data.
### only works when there are both left and right censored data.
    xvec <- as.vector(x)
    nn <- length(xvec)
    if (nn <= 1) 
        stop("Need more observations")
    if (length(d) != nn) 
        stop("length of x and d must agree")
    if (any((d != 0) & (d != 1) & (d != 2))) 
      stop("d must be 0(right-censored) or 1(uncensored) or 2(left-censored)")
    if (!is.numeric(xvec)) 
        stop("x must be numeric") 
    temp <- Wdataclean3(z=xvec, d=d, zc=zc, wt=wt)
    x <- temp$value
    d <- temp$dd
    w <- temp$weight
    INDEX10 <- which(d != 2)
    d[INDEX10[length(INDEX10)]] <- 1
    INDEX12 <- which(d != 0)
    d[INDEX12[1]] <- 1
    xd1 <- x[d == 1]
    if (length(xd1) <= 1) 
        stop("need more distinct uncensored obs.")
    xd0 <- x[d == 0]
    xd2 <- x[d == 2]
    wd1 <- w[d == 1]
    wd0 <- w[d == 0]
    wd2 <- w[d == 2]
    m <- length(xd0)
    mleft <- length(xd2)
    if ((m > 0) && (mleft > 0)) {
        pnew <- wd1/sum(wd1)
        n <- length(pnew)
        k <- rep(NA, m)
        for (i in 1:m) {
            k[i] <- 1 + n - sum(xd1 > xd0[i])
        }
        kk <- rep(NA, mleft)
        for (j in 1:mleft) {
            kk[j] <- sum(xd1 < xd2[j])
        }
        num <- 1
        while (num < maxit) {
            wd1new <- wd1
            sur <- cumsumsurv(pnew)   ## rev(cumsum(rev(pnew)))   3/2015 MZ
            cdf <- 1 - c(sur[-1], 0)
            for (i in 1:m) {
            wd1new[k[i]:n] <- wd1new[k[i]:n] + wd0[i] * pnew[k[i]:n]/sur[k[i]]
            }
            for (j in 1:mleft) {
        wd1new[1:kk[j]] <- wd1new[1:kk[j]] + wd2[j] * pnew[1:kk[j]]/cdf[kk[j]]
            }
            pnew <- wd1new/sum(wd1new)
            num <- num + 1
        }
        sur <- cumsumsurv(pnew)   ## rev(cumsum(rev(pnew)))   3/2015 MZ
        cdf <- 1 - c(sur[-1], 0)
 logel<-sum(wd1*log(pnew))+sum(wd0*log(sur[k])) + sum(wd2*log(cdf[kk]))
    }
    return(list(logEL=logel, time=xd1, jump=pnew, surv=1-cdf, prob=cdf))
}
### compute the self-consistent estimator without puting a mass at (0,2)
### Used in el.cen.EM/2 to get logel00, but also works standalone.
