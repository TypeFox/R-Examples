el.cen.EM.kmc<-function (x, d, fun = function(t) {

  t

}, mu, maxit = 25, error = 1e-09,debug.kmc=F, ...) 

{

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

  if (length(mu) != 1) 

    stop("check the dim of constraint mu")

  temp <- Wdataclean2(xvec, d)

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

  funxd1 <- fun(xd1, ...)

  xd0 <- x[d == 0]

  xd2 <- x[d == 2]

  wd1 <- w[d == 1]

  wd0 <- w[d == 0]

  wd2 <- w[d == 2]

  m <- length(xd0)

  mleft <- length(xd2)

  if ((m > 0) && (mleft > 0)) {

    pnew <- el.test.wt(funxd1, wt = wd1, mu)$prob

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

      sur <- rev(cumsum(rev(pnew)))

      cdf <- 1 - c(sur[-1], 0)

      for (i in 1:m) {

        wd1new[k[i]:n] <- wd1new[k[i]:n] + wd0[i] * pnew[k[i]:n]/sur[k[i]]

      }

      for (j in 1:mleft) {

        wd1new[1:kk[j]] <- wd1new[1:kk[j]] + wd2[j] * 

          pnew[1:kk[j]]/cdf[kk[j]]

      }

      temp8 <- el.test.wt(funxd1, wt = wd1new, mu)

      if (sum(abs(pnew-temp8$prob) )<error) break;

      pnew <- temp8$prob

      lam <- temp8$lam

      num <- num + 1

    }

    logel <- sum(wd1 * log(pnew)) + sum(wd0 * log(sur[k])) + 

      sum(wd2 * log(cdf[kk]))

    tempDB <- WCY(x = x, d = d, wt = w)

    logel00 <- tempDB$logEL

    funNPMLE <- sum(fun(tempDB$time) * tempDB$jump)

  }

  if ((m > 0) && (mleft == 0)) {

    temp3 <- WKM(x = x, d = d, w = w)

    logel00 <- temp3$logel

    funNPMLE <- sum(funxd1 * temp3$jump[temp3$jump > 0])

    pnew <- el.test.wt(x = funxd1, wt = wd1, mu = mu)$prob

    n <- length(pnew)

    k <- rep(NA, m)

    for (i in 1:m) {

      k[i] <- 1 + n - sum(xd1 > xd0[i])

    }

    num <- 1

    while (num < maxit) {

      wd1new <- wd1

      sur <- rev(cumsum(rev(pnew)))

      for (i in 1:m) {

        wd1new[k[i]:n] <- wd1new[k[i]:n] + wd0[i] * pnew[k[i]:n]/sur[k[i]]

      }

      temp9 <- el.test.wt(funxd1, wt = wd1new, mu)

      if (sum(abs(pnew-temp9$prob) )<error) break;

      pnew <- temp9$prob

      lam <- temp9$lam

      num <- num + 1

    }

    if (debug.kmc)cat('EM run time',num,'\n')

    sur <- rev(cumsum(rev(pnew)))

    logel <- sum(wd1 * log(pnew)) + sum(wd0 * log(sur[k]))

  }

  if ((m == 0) && (mleft > 0)) {

    kk <- rep(NA, mleft)

    for (j in 1:mleft) {

      kk[j] <- sum(xd1 < xd2[j])

    }

    pnew <- el.test.wt(funxd1, wt = wd1, mu)$prob

    n <- length(pnew)

    num <- 1

    while (num < maxit) {

      wd1new <- wd1

      cdf <- cumsum(pnew)

      for (j in 1:mleft) {

        wd1new[1:kk[j]] <- wd1new[1:kk[j]] + wd2[j] * 

          pnew[1:kk[j]]/cdf[kk[j]]

      }

      temp7 <- el.test.wt(funxd1, wt = wd1new, mu)

      if (sum(abs(pnew-temp7$prob) )<error) break;

      pnew <- temp7$prob

      lam <- temp7$lam

      num <- num + 1

    }

    logel <- sum(wd1 * log(pnew)) + sum(wd2 * log(cdf[kk]))

    dleft <- d

    dleft[dleft == 2] <- 0

    templeft <- WKM(x = -x, d = dleft, w = w)

    logel00 <- templeft$logel

    funNPMLE <- NA

  }

  if ((m == 0) && (mleft == 0)) {

    funNPMLE <- sum(funxd1 * wd1/sum(wd1))

    logel00 <- sum(wd1 * log(wd1/sum(wd1)))

    temp6 <- el.test.wt(funxd1, wt = wd1, mu)

    pnew <- temp6$prob

    lam <- temp6$lam

    logel <- sum(wd1 * log(pnew))

  }

  tval <- 2 * (logel00 - logel)

  list(loglik = logel, times = xd1, prob = pnew, funMLE = funNPMLE, 

       lam = lam, `-2LLR` = tval, Pval = 1 - pchisq(tval, df = 1))

}



## EL.CEN.EM2





el.cen.EM2.kmc<-function (x, d, xc = 1:length(x), fun, mu, maxit = 200, error = 1e-09,debug.kmc=F

         , ...) 

{

  xvec <- as.vector(x)

  d <- as.vector(d)

  mu <- as.vector(mu)

  xc <- as.vector(xc)

  n <- length(d)

  if (length(xvec) != n) 

    stop("length of d and x must agree")

  if (length(xc) != n) 

    stop("length of xc and d must agree")

  if (n <= 2 * length(mu) + 1) 

    stop("Need more observations")

  if (any((d != 0) & (d != 1) & (d != 2))) 

    stop("d must be 0(right-censored) or 1(uncensored) or 2(left-censored)")

  if (!is.numeric(xvec)) 

    stop("x must be numeric")

  if (!is.numeric(mu)) 

    stop("mu must be numeric")

  funx <- as.matrix(fun(xvec, ...))

  pp <- ncol(funx)

  if (length(mu) != pp) 

    stop("length of mu and ncol of fun(x) must agree")

  temp <- Wdataclean5(z = xvec, d, zc = xc, xmat = funx)

  x <- temp$value

  d <- temp$dd

  w <- temp$weight

  funx <- temp$xxmat

  INDEX10 <- which(d != 2)

  d[INDEX10[length(INDEX10)]] <- 1

  INDEX12 <- which(d != 0)

  d[INDEX12[1]] <- 1

  xd1 <- x[d == 1]

  if (length(xd1) <= 1) 

    stop("need more distinct uncensored obs.")

  funxd1 <- funx[d == 1, ]

  xd0 <- x[d == 0]

  xd2 <- x[d == 2]

  wd1 <- w[d == 1]

  wd0 <- w[d == 0]

  wd2 <- w[d == 2]

  m <- length(xd0)

  mleft <- length(xd2)

  if ((m > 0) && (mleft > 0)) {

    pnew <- el.test.wt2(x = funxd1, wt = wd1, mu = mu)$prob

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

      sur <- rev(cumsum(rev(pnew)))

      cdf <- 1 - c(sur[-1], 0)

      for (i in 1:m) {

        wd1new[k[i]:n] <- wd1new[k[i]:n] + wd0[i] * pnew[k[i]:n]/sur[k[i]]

      }

      for (j in 1:mleft) {

        wd1new[1:kk[j]] <- wd1new[1:kk[j]] + wd2[j] * 

          pnew[1:kk[j]]/cdf[kk[j]]

      }

      temp8 <- el.test.wt2(x = funxd1, wt = wd1new, mu = mu)

      pnew <- temp8$prob

      lam <- temp8$lambda

      num <- num + 1

    }

    logel <- sum(wd1 * log(pnew)) + sum(wd0 * log(sur[k])) + 

      sum(wd2 * log(cdf[kk]))

    logel00 <- WCY(x = x, d = d, wt = w)$logEL

  }

  ssst=1;

  if ((m > 0) && (mleft == 0)) {

    if (ssst==1){poo=0}else{poo=pnew} 

    pnew <- el.test.wt2(x = funxd1, wt = wd1, mu = mu)$prob

    n <- length(pnew)

    k <- rep(NA, m)

    for (i in 1:m) {

      k[i] <- 1 + n - sum(xd1 > xd0[i])

    }

    num <- 1

    while (num < maxit) {

      wd1new <- wd1

      sur <- rev(cumsum(rev(pnew)))

      for (i in 1:m) {

        wd1new[k[i]:n] <- wd1new[k[i]:n] + wd0[i] * pnew[k[i]:n]/sur[k[i]]

      }

      temp9 <- el.test.wt2(x = funxd1, wt = wd1new, mu = mu)

      pnew <- temp9$prob

      lam <- temp9$lambda

      num <- num + 1

    }

    sur <- rev(cumsum(rev(pnew)))

    logel <- sum(wd1 * log(pnew)) + sum(wd0 * log(sur[k]))

    logel00 <- WKM(x = x, d = d, zc = xc, w = w)$logel

    if ( sum(abs(poo-pnew))<error) break;

    ssst=ssst+1;

  }

  if (debug.kmc) {if (num == maxit) warning('EM may NOT converge!');}

  if ((m == 0) && (mleft > 0)) {

    kk <- rep(NA, mleft)

    for (j in 1:mleft) {

      kk[j] <- sum(xd1 < xd2[j])

    }

    pnew <- el.test.wt2(x = funxd1, wt = wd1, mu = mu)$prob

    n <- length(pnew)

    num <- 1

    while (num < maxit) {

      wd1new <- wd1

      cdf <- cumsum(pnew)

      for (j in 1:mleft) {

        wd1new[1:kk[j]] <- wd1new[1:kk[j]] + wd2[j] * 

          pnew[1:kk[j]]/cdf[kk[j]]

      }

      temp7 <- el.test.wt2(x = funxd1, wt = wd1new, mu = mu)

      pnew <- temp7$prob

      lam <- temp7$lambda

      num <- num + 1

    }

    logel <- sum(wd1 * log(pnew)) + sum(wd2 * log(cdf[kk]))

    dleft <- d

    dleft[dleft == 2] <- 0

    templeft <- WKM(x = -x, d = dleft, zc = xc, w = w)

    logel00 <- templeft$logel

  }

  if ((m == 0) && (mleft == 0)) {

    num <- 0

    temp6 <- el.test.wt2(x = funxd1, wt = wd1, mu)

    pnew <- temp6$prob

    lam <- temp6$lambda

    logel <- sum(wd1 * log(pnew))

    logel00 <- sum(wd1 * log(wd1/sum(wd1)))

  }

  tval <- 2 * (logel00 - logel)

  list(loglik = logel, times = xd1, prob = pnew, lam = lam, 

       iters = num, `-2LLR` = tval, Pval = 1 - pchisq(tval, 

                                                      df = length(mu)))

}





