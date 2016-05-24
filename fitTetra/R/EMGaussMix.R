EMGaussMix <-
function(y, ng=5, ptype="p.free", mutype=0, sdtype="sd.const", sd.fixed=0.05,
    p.start=rep(1/ng,ng), mu.start=seq(0.2, pi/2-0.2, by=(pi/2-0.4)/(ng-1)), sd.start=rep(0.075,ng),
    npar, maxiter, maxn.bin, nbin)  {
  
  #y is on the transformed (arcsine-square root) scale, like mu.start and sd.start

  isna <- is.na(y)
  yw <- y[!isna]
  n <- length(yw)
    binned <- F

  w <- rep(1,length(yw))

  if (n > maxn.bin) {
    binned <- T
    n.iv <- nbin
    lgt.iv <- (pi/2) / n.iv
    brks <- seq(0, pi/2, lgt.iv) 
    mids <- 0.5*lgt.iv + brks[-(n.iv+1)]
    brks[n.iv+1] <- brks[n.iv+1]+0.0001 # make sure that last break above pi/2
    ww <- hist(yw, brks, plot=F)$counts

    bins.filled  <- ww != 0
    n.iv2 <- sum(bins.filled)

    mis.bin <-(1:n.iv)[!bins.filled]
    binnrs <- numeric(n.iv)
    binnrs[mis.bin] <- NA
    binnrs[!is.na(binnrs)] <- 1:n.iv2

    #yw.mmbr <- floor(yw/lgt.iv)+1 #to which bin (in 1:nbin) does each yw belong? 
                                   #wrong: assumes hist intervals are closed at bottom, but they are closed at top
    yw.mmbr <- ceiling(yw/lgt.iv)  #to which bin (in 1:nbin) does each yw belong?
    yw.mmbr[yw.mmbr==0] <- 1       #else if yw==0 a non-existing bin 0 would be selected
    yw.mmbr <- binnrs[yw.mmbr]     #to which bin (in 1:n.iv2) does each yw belong?

    yw <- mids[bins.filled]    # empty bins are removed
    w  <- ww[bins.filled]      # weight (number of samples) in each bin
  }
  W <- diag(w)

  result <- list (loglik=NA, npar=npar, AIC=NA, BIC=NA, psi=NA, post=NA, nobs=n, iter=NA, message="")
  psi <- list(mu=mu.start, sigma=sd.start, p=p.start/sum(p.start))
  o <- order(psi$mu); psi <- orderpsi(psi,o)
  llnw <- loglikf(yw,w,psi) #llnw is the new total loglik of y with frequencies w given the mu, sigma and p vectors in psi
  iter <- 0
  tryCatch( {
    repeat {
      iter <- iter+1
      ll <- llnw
      Z <- EMGaussExp.vectorized(yw, psi) #per yw ng numbers: probs that y belongs to each peak

      result$message <- tryCatch( {
        psi <- EMGaussMax(Z, W, yw, ptype, mutype, sdtype, sd.fixed, p.start)

        ""
      }, error=function(ex) { paste(ex$message,"in EMGaussMix") } )
      if (result$message!="") break
        #Note: on error here control jumps back to CodomMarker, not caught by local tryCatch !?
      llnw <- loglikf(yw,w,psi)
      if (is.na(llnw) || (max(llnw - ll) < 0.0000001) ||
           (maxiter>0 && iter>maxiter) ) break
    }
    if (result$message=="") {
      if ( is.na(llnw) ) result$message<-"llnw==NA in EMGaussMix"
      else if (max(llnw - ll) < 0.0000001) {
        o <- order(psi$mu); psi <- orderpsi(psi,o); Z <- Z[,o]

        NAmat <- matrix(rep(!isna,ng), ncol=ng,  byrow=F)
        z <- matrix(nrow=length(y), ncol=ng)
        if (binned) {
          z[NAmat] <- Z[yw.mmbr,]
        }
        else z[NAmat] <- Z

        result$post <- z
        llnw <- loglikf(y[!isna],rep(1,length(y[!isna])),psi)
        result$loglik <- llnw
        result$AIC <- -2*llnw + 2*npar
        result$BIC <- -2*llnw + log(n)*npar
        result$psi <- psi
        result$iter <- iter
      }
      else result$message <- "iter>maxiter in EMGaussMix"
    }
  }, error=function(ex) { result$message <- paste(ex$message," in EMGaussMix",sep="") } )
  result
}
