lgarchSim <-
function(n, constant=0, arch=0.05, garch=0.9, xreg=NULL,
  backcast.values=list(lnsigma2=NULL, lnz2=NULL, xreg=NULL),
  check.stability=TRUE, innovations=NULL, verbose=FALSE,
  c.code=TRUE)
{
  #check arguments:
  if(is.null(constant)) constant <- 0
  if(is.null(arch)) arch <- 0
  if(is.null(garch)) garch <- 0
  #if(is.null(asym)) asym <- 0

  #orders:
  arch.p <- length(arch)
  garch.q <- length(garch)
  maxpq <- max(arch.p, garch.q)
  nmaxpq <- n+maxpq
  #### NOTE: asym not implemented yet
  #asym.r <- length(asym)
  #maxpqr <- max(arch.p, garch.q, asym.r)
  #nmaxpqr <- n+maxpqr
  ####
  #### NOTE: c.code not implemented yet

  #make phi:
  phi <- c(arch,rep(0,maxpq-arch.p)) + c(garch,rep(0,maxpq-garch.q))

  #check stability:
  if(check.stability){
    roots <- polyroot(c(1,-phi))
    if(any(abs(roots) <= 1)){
      mssg <- paste("NOTE: The log-garch model may not be stable (one or more AR-roots is on or inside the unit circle)")
      print(mssg)
    }
  }

  #xreg, zoo-index:
  if(is.null(xreg)){
    result.index <- 1:n
    xreg <- rep(0,nmaxpq)
    xregMean <- 0
  }else{
    xreg <- as.zoo(xreg)
    result.index <- index(xreg)
    xreg <- coredata(xreg)
    if(is.null(backcast.values$xreg)){
      xregMean <- mean(xreg)
      xreg.init <- rep(xregMean,maxpq)
    }else{
      xreg.init <- backcast.values$xreg
    }
    xreg <- c(xreg.init, xreg)
  }

  #innovation:
  if(is.null(innovations)){
    z <- rnorm(n)
  }else{
    z <- innovations
  }
  z2 <- z^2
  lnz2 <- log(z2)
  if(is.null(backcast.values$lnz2)){
    Elnz2 <- mean(lnz2)
    lnz2 <- c(rep(Elnz2, maxpq), lnz2)
  }else{
    lnz2 <- c(backcast.values$lnz2, lnz2)
  }
  mLnz2 <- NULL
  for(k in 1:arch.p){
    mLnz2 <- cbind(mLnz2, glag(lnz2, k=k, pad=TRUE,
      pad.value=Elnz2))
  }

  #innov series:
  if(arch.p > 0){
    lnz2.innov <- as.vector(mLnz2 %*% arch)
  }else{ lnz2.innov <- rep(0,nmaxpq) }
  innov <- constant + lnz2.innov + xreg

  #lnsigma2:
  if(is.null(backcast.values$lnsigma2)){
    Elnsigma2 <- (constant+sum(arch)*Elnz2+xregMean)/(1-sum(phi))
    if(abs(Elnsigma2)==Inf){
      mssg <- paste("NOTE: The backcast value(s) of lnsigma2 is Inf")
      print(mssg)
    }
    lnsigma2 <- c(rep(Elnsigma2,maxpq),rep(0,n))
  }else{
    lnsigma2 <- c(backcast.values$lnsigma2,rep(0,n))
  }

  #recursion:
  maxpq1 <- maxpq+1
  phisum <- rep(0,nmaxpq)
  if(c.code){
    tmp <- LGARCHSIM(maxpq, nmaxpq, lnsigma2, phi, phisum, innov)
    lnsigma2 <- tmp$lnsigma2
  }else{
    for(i in maxpq1:nmaxpq){
      phisum[i] <- sum(phi*lnsigma2[c(i-1):c(i-maxpq)])
      lnsigma2[i] <- innov[i] + phisum[i]
    } #end for loop
  } #end if(c.code)

  #output:
  lnsigma2 <- lnsigma2[-c(1:maxpq)] #rm initial vals
  if(verbose){
    sigma <- exp(lnsigma2/2)
    y <- sigma*z
    lny2 <- log(y^2)
    lnz2 <- log(z^2)
    result <- cbind(y, lny2, sigma, lnsigma2, z, lnz2)
    colnames(result)[3:4] <- c("sd","lnsd2")
  }else{
    sigma <- exp(lnsigma2/2)
    result <- sigma*z
  }
  result <- zoo(result, order.by=result.index)
  return(result)
}
