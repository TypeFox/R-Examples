weights.negbin <- function(y, mu, alpha=1, bw=NULL, raf='SCHI2', tau=NULL, nrep=100, type=c('unif', 'norm')) {
  pesi <- function(x, raf, tau=0.1) {
    # x: residui di Pearson
    gkl <- function(x, tau) {
      if (tau!=0)
        x <- log(tau*x+1)/tau
      return(x)
    }
    pwd <- function(x, tau) {
      if(tau==Inf)
        x <- log(x+1)
      else
        x <- tau*((x + 1)^(1/tau) - 1)
      return(x)
    }
    x <- ifelse(x < 1e-10, 0, x)
    ww <- switch(raf,
       ##park+basu+2003.pdf
       GKL = gkl(x, tau),
       ##lindsay+1994.pdf                 
       PWD = pwd(x, tau),
       HD =  2*(sqrt(x + 1) - 1) ,
       NED =  2 - (2 + x)*exp(-x) ,
       SCHI2 =  1-(x^2/(x^2 +2)) )

    if (raf!='SCHI2') {
      ww <- (ww + 1)/(x + 1)
    }
    ww[ww > 1] <- 1
    ww[ww < 0] <- 0
    ww[is.infinite(x)] <- 0
    return(ww)
  }
  type <- match.arg(type)
  n <- length(y)
  size <- 1/alpha
  p <- q <- rep(0, n)
  qmat <- mmat <- matrix(0, n, nrep)
  for (i in 1:nrep) {
    u <- runif(n)
    for (k in 1:n) {
      p[k] <- pnbinom(y[k], size=size, mu=mu[k])
      v <- ifelse(y[k] <= size, dnbinom(y[k], size=size, mu=mu[k])*u[k], 0)
      q[k]  <- p[k]-v
    }
    if (type=='unif') {
#    qecdf <- ecdf(q)
      qecdf <- approxfun(x=sort(q), y=((1:n)-0.5)/n, rule=2)
      qmat[,i] <- qecdf(q+bw)-qecdf(q-bw)
      mmat[,i] <- punif(q+bw)-punif(q-bw)
    } else {
      q <- qnorm(q)
      qden <- density(q, bw=bw)
      qden <- approxfun(x=qden$x, y=qden$y, rule=2)
      qmat[,i] <- qden(q)
      mmat[,i] <- dnorm(q, mean=0, sd=(1+sqrt(bw)))
    }
  }
  delta <- qmat/mmat-1
  delta <- apply(delta, 1, mean)
  w <- pesi(x=delta, raf=raf, tau=tau)
  return(w)
}
