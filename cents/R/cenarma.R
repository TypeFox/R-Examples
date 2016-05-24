#Source: cenarma.R
#censored arma fitting function
#
cenarma <- function(y, iy, p=0, q=0, include.mean=TRUE, 
                    verbose=FALSE, MaxIter=100, ETOL=1e-5,
                    algorithm=c("exact","approx"), ...){
#if p=q=0, ignore iy for the time being. Later include censoring!!!
#
  if (missing(iy))
    iy <- ifelse(is.na(y), "na", "o")
  stopifnot(length(iy)==length(y))
  alg <- match.arg(algorithm)
  n <- length(y)
  muk2 <- numeric(n)
  yr <- rev(y)
  iyr <- rev(iy)
  nobs <- sum(as.integer(iy=="o"))
  muHat <- ifelse(include.mean, mean(y[iy=="o"]), 0) 
  sigysqHat <- sum((y[iy=="o"]-muHat)^2)/nobs
  sigyHat <- sqrt(sigysqHat)
  #initial estimates obtained by treating censored values correct values
  ans <- try(arima(y, order=c(p,0,q), include.mean=include.mean, ...), silent=TRUE)
  if (identical(class(ans), "try-error")) {
    warning("error occurred in arima")
    return(NA)
  }
  if (ans$code!=0) 
    warning("arima returned convergence code: ", ans$code)
  betaHat0 <- coef(ans)
  arHat <- maHat <- numeric(0)
  if(p!=0) arHat <- betaHat0[1:p]
  if(q!=0) maHat <-  betaHat0[(p+1):(p+q)]
  muHat <- ifelse(include.mean, betaHat0[p+q+1], 0)
  sigmaSqAHat <- ans$sigma2
  logL0 <- -Inf
  err <- 1e3
  iter <- 0
  while(abs(err)>ETOL && iter<MaxIter) {
    iter <- iter+1
    r <- switch(alg,
      approx = sigysqHat*ARMAacf(ar = arHat, ma = maHat, lag.max = n-1),
      exact = tacvfARMA(phi=arHat, theta=-maHat, maxlag=n-1, sigma2=sigmaSqAHat)
      )
    x <- y
    muk2[1] <- muHat
    x[1] <- switch(iy[1], o=y[1], na=muHat, L=EZR(muHat,sigyHat,y[1]), 
                   R=EZL(muHat,sigyHat,y[1]))
    phi <- r[2]/r[1]
    sigmasqk <- r[1] * (1 - phi^2)
    muk <- muHat + phi * (x[1]-muHat)
    muk2[2] <- muk
    x[2] <- switch(iy[2], o=y[2], na=muk, L=EZR(muk,sqrt(sigmasqk),y[2]), 
                   R=EZL(muk,sqrt(sigmasqk),y[2]))
    sigmasqkm1 <- sigmasqk
    for (k in 2:(n - 1)) {
      if (sigmasqkm1 < 0) 
      	sigmasqkm1 <- 0
      phikk <- (r[k + 1] - phi %*% rev(r[2:k]))/sigmasqkm1
      sigmasqk <- sigmasqkm1 * (1 - phikk^2)
      phinew <- phi - phikk * rev(phi)
      phi <- c(phinew, phikk)
      sigmasqkm1 <- sigmasqk
      muk <- muHat+crossprod(phi, rev(x[1:k]-muHat))
      muk2[k+1] <- muk
      x[k+1] <- switch(iy[k+1], o=y[k+1], na=muk, L=EZR(muk,sqrt(sigmasqk),y[k+1]), 
                       R=EZL(muk,sqrt(sigmasqk),y[k+1]))
    }
    ans <- try(arima(x, order=c(p,0,q), include.mean=include.mean, ...), silent=TRUE)
    if (identical(class(ans), "try-error")) 
      return(NA)
    if (ans$code!=0)
      warning("arima returned convergence code: ", ans$code)
    sigmaSqAHat <- ans$sigma2 
    logL <- ans$loglik
    betaHat <- coef(ans)
    arHat <- maHat <- numeric(0)
    if(p!=0) arHat <- betaHat[1:p]
    if(q!=0) maHat <-  betaHat[(p+1):(p+q)]
    muHat <- ifelse(include.mean, betaHat[p+q+1], 0)
    sigysqHat <- sum((y[iy=="o"]-muHat)^2)/nobs
    sigyHat <- sqrt(sigysqHat)
    err <- logL - logL0
    logL0 <- logL
    #
    if (verbose) {
      cat("\nIteration: ", iter, fill=TRUE)
      cat(paste0("logL = ", logL), fill=TRUE)
      cat(paste0("err = ", err), fill=TRUE)      
      cat(paste0("muHat = ", muHat), fill=TRUE)
      cat(paste0("arHat = ", arHat), fill=TRUE)
      cat(paste0("maHat = ", maHat), fill=TRUE)
      cat(paste0("sigmasqk = ", sigmasqk), fill=TRUE)
      }
    }
  se1 <- sqrt(sigmaSqAHat/(var(x)*n))
  out <- list(outarima=ans, p=p, q=q, include.mean=include.mean, se1=se1, 
       dataSummary=t(as.matrix(table(iy))),
       exitStatus=ifelse(err<ETOL,"converged","MaxIter reached"),
       y=y, iy=iy, call = match.call())
  class(out) <- "Cenarma"
  out
}

