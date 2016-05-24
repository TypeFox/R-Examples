##charles
## Some minor changes to functions in the logcondens package.


## This is modified from original code to allow for 'phi' to be as in the
## laplace distribution, e.g. arbitrarily large knot lengths.  We do not
## yet accomodate either of x or y here being large positive (where "x
## large" means exp(x)=Inf).  This is unexpected for us generally, and
## would require |y-x| to be very small to allow the density integral to be
## 1.  Thus it cannot appropriately be handled here, but would require
## processing in Local_LL_all.  Alternatively data could just be rescaled.
## THUS we assume: exp(x) and exp(y) are computable, potentially very small
## or 0.  exp|y-x| may be infinite.  Also relevant: it is not quite true
## that exp(r)==Inf implies exp(-r)==0. e.g. r=721
## Also: I'm not sure what the point of the middle case is.  Left it in
## because that was how it was initially computed, not sure if it's better
## or not.
## The modified J10,J20,J11 give identical answers to the old versions when
## the old versions don't return Inf or NaN.


## J10 <- function (x, y) 
## {
##     m <- length(x)
##     z <- exp(x)
##     d <- y - x
##     LARGE <- (abs(d) > 0.01)
##     SMALL <- !LARGE
##     II <- (1:m)[LARGE]
##     z[II] <-  (exp(y[II]) - z[II] - d[II]*z[II])/(d[II]^2)    
##     II <- (1:m)[SMALL]
##     z[II] <- z[II] * (1/2 + d[II] * (1/6 + d[II] * (1/24 + d[II] * 
##         (1/120 + d[II]/720))))#by defn II, z[II] hasn't been modified yet
##     return(z)
## }



## J20 <- function (x, y) 
## {
##     m <- length(x)
##     z <- exp(x)
##     d <- y - x
##     ##ed <- exp(d)
##     ##LARGE <- ed==Inf
##     ##MED <- (abs(d) > 0.02) & !LARGE
##     ##SMALL <- !(MED | LARGE)
##     LARGE <- (abs(d)>0.02)
##     SMALL <- !LARGE
##     II <- (1:m)[LARGE]
##     z[II] <- 2* (exp(y[II]) - z[II] - z[II]*d[II] - z[II]*d[II]^2/2 ) / (d[II]^3)
##     ##II <- (1:m)[MED]
##     ##z[II] <- 2 * z[II] * (exp(d[II]) - 1 - d[II] - d[II]^2/2)/(d[II]^3)
##     II <- (1:m)[SMALL]
##     z[II] <- z[II] * (1/3 + d[II] * (1/12 + d[II] * (1/60 + d[II] * 
##         (1/360 + d[II]/2520))))
##     return(z)
## }

## J11 <-   function (x, y) 
## {
##   m <- length(x)
##   z <- exp(x)
##   d <- y - x
##   ##ed <- exp(d)
##   ##LARGE <- ed==Inf
##   ##MED <- (abs(d) > 0.02) & !LARGE
##   ##SMALL <- !(MED | LARGE)
##   LARGE <- (abs(d)>.02)
##   SMALL <- !LARGE
##   II <- (1:m)[LARGE]
##   z[II] <- (d[II]*(exp(y[II])+z[II]) - 2*(exp(y[II])-z[II])) / (d[II]^3)
##   ##II <- (1:m)[MED]
##   ##z[II] <- z[II] * (d[II] * (exp(d[II]) + 1) -
##   ##                  2 * (exp(d[II]) -  1))/(d[II]^3)
##   II <- (1:m)[SMALL]
##   z[II] <- z[II] * (1/6 + d[II] * (1/12 + d[II] * (1/40 + d[II] * 
##                                                    (1/180 + d[II]/1008))))
##   return(z)
## }



## J00 <- function (x, y, v = 1) 
## {
##     m <- length(x)
##     z <- exp(x)
##     d <- y - x
##     ed <- exp(v*d)
##     LARGE <- (abs(d)>0.005)
##     SMALL <- !LARGE
##     II <- (1:m)[LARGE]
##     z[II] <- (exp(v*y[II]+(1-v)*x[II]) - z[II]) / d[II]
##     II <- (1:m)[SMALL]
##     z[II] <- z[II] * (v + d[II] * (v/2 + d[II] * (v/6 + d[II] * 
##         (v/24 + d[II] * v/120))))
##     return(z)
## }



## Only a very slight modification from the original.  Added the
## "if(length(LIs)>0)" bit for error checking.  In the future one could
## accomodate exp(phi)=Inf because necessarily dx will be tiny soince phi
## concave and int exp(phi)=1.  Seems unlikely to be useful at this point,
## though.
Local_LL_all <- function (x, w, phi) 
{
  n <- length(x)
  dx <- diff(x)
  ll <- sum(w * phi) - sum(dx * J00(phi[1:(n - 1)], phi[2:n]))
  LIs <- (1:n)[exp(phi)==Inf] ##Large Indices
  if (length(LIs) > 0 ) stop("Local_LL_all Error: We have not yet accounted for extraordinarily steep/large phi.  This is probably an error.")
  grad <- matrix(w, ncol = 1)
  grad[1:(n - 1)] <- grad[1:(n - 1)] - (dx * J10(phi[1:(n -  1)], phi[2:n]))
  grad[2:n] <- grad[2:n] - (dx * J10(phi[2:n], phi[1:(n - 1)]))  
  tmp <- c(dx * J20(phi[1:(n - 1)], phi[2:n]), 0) +
    c(0, dx *  J20(phi[2:n], phi[1:(n - 1)]))
  tmp <- tmp + mean(tmp) * 1e-12
  mhess2 <- matrix(0, nrow = n, ncol = n)
  mhess3 <- mhess2
  mhess1 <- tmp
  tmp <- c(0, dx * J11(phi[1:(n - 1)], phi[2:n]))
  tmp.up <- diag(tmp[2:n], nrow = n - 1, ncol = n - 1)
  mhess2[1:(n - 1), 2:n] <- tmp.up
  mhess3[2:n, 1:(n - 1)] <- diag(tmp[2:n], nrow = n - 1, ncol = n -  1)
  mhess <- diag(mhess1) + mhess2 + mhess3
  phi_new <- phi + solve(mhess) %*% grad
  dirderiv <- t(grad) %*% (phi_new - phi)
  return(list(ll = ll, phi_new = phi_new, dirderiv = dirderiv))
}

## why isn't this named 'LocalExtend.mode"

## as usual, 'constr' is 2 consecutive indices in x2, the knots,
## corresponding to the two constrained ones; or, if no constraint,
## it can be c(i,i) for any i.
LocalExtend <- function (x, IsKnot, x2, phi2, constr=NULL) {
    n <- length(x)
    K <- (1:n) * IsKnot
    K <- K[K > 0]
    phi <- 1:n * 0
    phi[K] <- phi2
    for (k in 1:(length(K) - 1)) {
        if (K[k + 1] > (K[k] + 1)) {
            ind <- (K[k] + 1):(K[k + 1] - 1)
            lambda <- (x[ind] - x2[k])/(x2[k + 1] - x2[k])
            phi[ind] <- (1 - lambda) * phi2[k] + lambda * phi2[k + 
                1]
        }
      }
    if (length(constr)>1 && constr[1] != constr[2]){
      ind <- K[constr[1]]:K[constr[2]] ##constr1,2 could be consecutive in x
      phi[ind] <- rep(phi2[constr[1]], length=length(ind))
    }
    return(matrix(phi, ncol = 1))
}


## prec defines when a divisor is assumed to be ==0.
intF <- function (s, x, phi, Fhat,
                  prec=1e-10 ) 
{
  ##x assumed to be sorted increasing.
  n <- length(x)
  lows <- s<x[1]
  upps <- s>x[n]
  dx <- c(NA, diff(x))
  dphi <- c(NA, diff(phi))
  f <- exp(phi)
  F <- Fhat
  intF.xi <- c(0, rep(NA, n - 1))
  for (i in 2:n) {
    if (abs(dphi[i]) < prec){ ## "==0" often fails.
    intF.xi[i] <- dx[i] * (F[i-1] + f[i-1]*dx[i]/2)
    }
    else{
      intF.xi[i] <- dx[i] * (F[i - 1] + dx[i]/dphi[i] * 
                             (J00(phi[i - 1], phi[i], 1) - f[i - 1]))
    }
  }
  intF.xi <- cumsum(intF.xi)
  intF.s <- rep(0, length(s))
  for (k in 1:length(s)) {
    if (!lows[k]){ ##leave at 0.
      extra <- 0
      if (!upps[k]) mys <- s[k]
      else {
        extra <- s[k]-x[n]
        mys <- x[n]
      }
      j <- max((1:n)[x <= mys])
      j <- min(j, n - 1)
      Fj <- F[j]
      ds <- (mys-x[j])
      if (abs(dphi[j+1]) < prec){
        intF.s[k] <- intF.xi[j] +  ds * (Fj  + ds*f[j+1]/2) 
      }
      else{
        intF.s[k] <- intF.xi[j] + ds * Fj + dx[j+1] *
          (dx[j + 1]/dphi[j + 1] *
           J00(phi[j], phi[j + 1], ds/(dx[j + 1])) -
           ds/(dphi[j + 1]) * f[j])
      }
      intF.s[k] <- intF.s[k]+extra
    }
  }
  ##intF.s[lows] <- rep(0,sum(lows))
  ##intF.s[upps] <- intF.s[upps]+s[upps]-x[n]
  return(intF.s)
}

## intECDF <- function (s, x) 
## {
##   n <- length(x)
##   lows <- s<x[1]
##   upps <- s>x[n]
##   n <- length(x)
##   x <- sort(x)
##   dx <- c(0, diff(x))
##   intED.xi <- cumsum(dx * ((0:(n - 1))/n))
##   intED.s <- rep(0, length(s))
##   for (k in 1:length(s)) {
##     if ( !lows[k]){
##       mys <- s[k]
##       if (upps[k]) mys <- x[n]
##       j <- max((1:n)[x <= mys])
##       j <- min(j, n - 1)
##       xj <- x[j]
##       intED.s[k] <- intED.xi[j] + (mys - xj) * j/n
##     }
##   }
##   ##intF.s[lows] <- rep(0,sum(lows))
##   intED.s[upps] <- intED.s[upps]+s[upps]-x[n]
##   return(intED.s)
## }



## activeSetLogCon <- function (x, w = NA, prec=10^-10, print = FALSE, logfile=NULL) 
## {
##   if (!is.null(logfile) && !is.na(logfile)) sink(NULL)
##   n <- length(x)
##   if (sum(x[2:n] <= x[1:(n - 1)]) > 0) {
##     cat("We need strictly increasing numbers x(i)!\n")
##   }
##   if (max(is.na(w)) == 1) {
##     w <- rep(1/n, n)
##   }
##   if (sum(w <= 0) > 0) {
##     cat("We need strictly positive weights w(i)!\n")
##   }
##   w <- w/sum(w)
##   phi <- LocalNormalize(x, 1:n * 0)
##   IsKnot <- 1:n * 0
##   IsKnot[c(1, n)] <- 1
##   res1 <- LocalMLE(x, w, IsKnot, phi, prec)
##   phi <- res1$phi
##   L <- res1$L
##   conv <- res1$conv
##   H <- res1$H
##   iter1 <- 1
##   while ((iter1 < 500) & (max(H) > prec * mean(abs(H)))) {
##     IsKnot_old <- IsKnot
##     iter1 <- iter1 + 1
##     tmp <- max(H)
##     k <- (1:n) * (H == tmp)
##     k <- min(k[k > 0])
##     IsKnot[k] <- 1
##     res2 <- LocalMLE(x, w, IsKnot, phi, prec)
##     phi_new <- res2$phi
##     L <- res2$L
##     conv_new <- res2$conv
##     H <- res2$H
##     while ((max(conv_new) > prec * max(abs(conv_new)))) {
##       JJ <- (1:n) * (conv_new > 0)
##       JJ <- JJ[JJ > 0]
##       tmp <- conv[JJ]/(conv[JJ] - conv_new[JJ])
##       lambda <- min(tmp)
##       KK <- (1:length(JJ)) * (tmp == lambda)
##       KK <- KK[KK > 0]
##       IsKnot[JJ[KK]] <- 0
##       phi <- (1 - lambda) * phi + lambda * phi_new
##       conv <- pmin(c(LocalConvexity(x, phi), 0))
##       res3 <- LocalMLE(x, w, IsKnot, phi, prec)
##       phi_new <- res3$phi
##       L <- res3$L
##       conv_new <- res3$conv
##       H <- res3$H
##       if (print == TRUE) {
##         print(paste("iter1=", iter1, " / L=", round(L, 
##                                                     4), " / max(H)=", round(max(H), 4), sep = ""))
##       }
##     }
##     phi <- phi_new
##     conv <- conv_new
##     if (sum(IsKnot != IsKnot_old) == 0) {
##       break
##     }
##     if (print == TRUE) {
##       print(paste("iter1=", iter1, " / L=", round(L, 4), 
##                   " / max(H)=", round(max(H), 4), sep = ""))
##     }
##   }
##   Fhat <- LocalF(x, phi)
##   phia <- max(phi)
##   aidx <- which(phi==phia)
##   aval <- x[aidx]
##   ## some of these return vals are to be compatible with the logcondens package.
##   res1 <- list(xn=x, ##CHANGE THIS 
##                x = x, w=w,
##                IsKnot = IsKnot,
##                L = L,                              
##                knots=x[IsKnot==1], ## this has been changed! in my old code it returned INDICES. now it returns VALUES
##                n=length(x),
##                ##m=, ## 'm' is length(x), n=length(unique(x))
##                phi = as.vector(phi), 
##                fhat=as.vector(exp(phi)),
##                Fhat = as.vector(Fhat),
##                H = as.vector(H),
##                a=list(val=aval,idx=aidx, isx=TRUE))
##   phi.f <- function(x0){
##     ##myf <- function(x00){evaluateLogConDens(x00,x=x,phi,Fhat,IsKnot)[1]}
##     #myf <- function(x00){evaluateLogConDens(xs,res=res1,which=1)}
##     #apply(matrix(x0),1,myf)
##     evaluateLogConDens(x0,res=res1,which=1)[,2] ##[,"log-density"]
##   }
##   fhat.f <- function(x0){
##     ##myf <- function(x00){evaluateLogConDens(x00,x=x,phi,Fhat, IsKnot)[2]}
##     ##apply(matrix(x0),1,myf)
##     evaluateLogConDens(x0,res=res1,which=2)[,3] 
##   }
##   Fhat.f <- function(x0){
##     ##myf <- function(x00){evaluateLogConDens(x00,x=x,phi,Fhat, IsKnot)[3]}
##     ##apply(matrix(x0),1,myf)
##     evaluateLogConDens(x0,res=res1,which=3)[,4]
##   }
##   E.f <- intFfn(x,phi,Fhat)
##   {##get phiP=phi prime= deriv of phi; left and right derivs are L and R.
##     KK <- (1:n)[as.logical(IsKnot)]
##     phiK <- phi[KK]
##     slopes <- diff(phiK)/diff(x[KK])
##     phiPR.f <- stepfun(x=x[KK], c(slopes[1],slopes,tail(slopes,1)), right=FALSE,f=0)
##     phiPL.f <- stepfun(x=x[KK], c(slopes[1],slopes,tail(slopes,1)), right=TRUE,f=1)
##     phiPL <- phiPL.f(x)
##     phiPR <- phiPR.f(x)
##   }
##   if (!is.null(logfile) && !is.na(logfile)) sink(NULL)
##   return(c(res1,
##            phi.f=phi.f, fhat.f=fhat.f, Fhat.f=Fhat.f, E.f=E.f,
##            phiPL=phiPL,phiPR=phiPR,
##            phiPL.f=phiPL.f,phiPR.f=phiPR.f))
## }


## Need to write wrapper function, logConDens? or does the
## function shipping with the logcondens package do the job still?
activeSetLogCon <-
  function(x, xgrid=NULL,print = FALSE,   w = NA, prec=10^-10) 
{
  if (print){ print("ASLC: Beginning")}
  ## if (!is.null(logfile) && !is.na(logfile)) sink(logfile)## just use sink() outside
  xn <- sort(x)
  if ((!identical(xgrid, NULL) & (!identical(w, NA)))) {
    stop("If w != NA then xgrid must be NULL!\n")
  }
  if (identical(w,NA)){
    tmp <- preProcess(x,xgrid=xgrid)
    x <- tmp$x
    w <- tmp$w
    sig <- tmp$sig ##nonpara est of sd.
  }
  else { ## if (!identical(w, NA)) {
    if (abs(sum(w) - 1) > prec) stop("activeSetLogCon Error: weights w do not sum to 1.")
    tmp <- cbind(x, w)
    tmp <- tmp[order(x), ]
    x <- tmp[, 1]
    w <- tmp[, 2]
    est.m <- sum(w * x)
    est.sd <- sum(w * (x - est.m)^2)
    est.sd <- sqrt(est.sd * length(x)/(length(x) - 1))
    sig <- est.sd
  }  
  n <- length(x)
  phi <- LocalNormalize(x, 1:n * 0)
  IsKnot <- 1:n * 0
  IsKnot[c(1, n)] <- 1
  res1 <- LocalMLE(x, w, IsKnot, phi, prec)
  phi <- res1$phi
  L <- res1$L
  conv <- res1$conv
  H <- res1$H
  iter1 <- 1
  while ((iter1 < 500) & (max(H) > prec * mean(abs(H)))) {
    IsKnot_old <- IsKnot
    iter1 <- iter1 + 1
    tmp <- max(H)
    k <- (1:n) * (H == tmp)
    k <- min(k[k > 0])
    IsKnot[k] <- 1
    res2 <- LocalMLE(x, w, IsKnot, phi, prec)
    phi_new <- res2$phi
    L <- res2$L
    conv_new <- res2$conv
    H <- res2$H
    while ((max(conv_new) > prec * max(abs(conv_new)))) {
      JJ <- (1:n) * (conv_new > 0)
      JJ <- JJ[JJ > 0]
      tmp <- conv[JJ]/(conv[JJ] - conv_new[JJ])
      lambda <- min(tmp)
      KK <- (1:length(JJ)) * (tmp == lambda)
      KK <- KK[KK > 0]
      IsKnot[JJ[KK]] <- 0
      phi <- (1 - lambda) * phi + lambda * phi_new
      conv <- pmin(c(LocalConvexity(x, phi), 0))
      res3 <- LocalMLE(x, w, IsKnot, phi, prec)
      phi_new <- res3$phi
      L <- res3$L
      conv_new <- res3$conv
      H <- res3$H
      if (print == TRUE) {
        print(paste("iter1=", iter1, " / L=", round(L, 4),
                    " / max(H)=", round(max(H), 4),
                    " / #knots = ", sum(IsKnot),
                    sep = ""))
      }
    }
    phi <- phi_new
    conv <- conv_new
    if (sum(IsKnot != IsKnot_old) == 0) {
      break
    }
    if (print == TRUE) {
      print(paste("iter1=", iter1, " / L=", round(L, 4), 
                  " / max(H)=", round(max(H), 4),
                  " / #knots = ", sum(IsKnot),
                  sep = ""))
    }
  }
  ## Now put work together for return values
  KK <- (1:n)[as.logical(IsKnot)]
  Fhat <- LocalF(x, phi)
  phia <- max(phi)
  aidx <- which(phi==phia)
  aval <- x[aidx]
  dlcMode <- list(val=aval,idx=aidx, isx=TRUE)
  class(dlcMode) <- "dlc.mode"
  ## some of these return vals are to be compatible with the logcondens package.
  res1 <- list(xn=xn, ## original passed-in x's, sorted. May contain repeats. 
               x=x, ##duplicates removed
               w=w,
               L = L,
               IsKnot = IsKnot,
               knots=x[IsKnot==1], ## this has been changed! in my old code it returned INDICES. now it returns VALUES
               phi = as.vector(phi), 
               fhat=as.vector(exp(phi)),
               Fhat = as.vector(Fhat),
               H = as.vector(H),
               ##a=list(val=aval,idx=aidx, isx=TRUE),
               n=length(xn),
               m=n,## m =length(x) <= length(xn)
               mode=aval, ## redundant, for backwards compatibility.
               ## note: "dlcMode" replaces old "a".
               ## dlcMode$val == mode
               dlcMode=dlcMode, ##naming: "mode" taken/exists already. "Mode" i'm against case-sens differences. 
               sig=sig)
  phi.f <- function(x0){
    ##myf <- function(x00){evaluateLogConDens(x00,x=x,phi,Fhat,IsKnot)[1]}
    #myf <- function(x00){evaluateLogConDens(xs,res=res1,which=1)}
    #apply(matrix(x0),1,myf)
    evaluateLogConDens(x0,res=res1,which=1)[,2] ##[,"log-density"]
  }
  fhat.f <- function(x0){
    ##myf <- function(x00){evaluateLogConDens(x00,x=x,phi,Fhat, IsKnot)[2]}
    ##apply(matrix(x0),1,myf)
    evaluateLogConDens(x0,res=res1,which=2)[,3] 
  }
  Fhat.f <- function(x0){
    ##myf <- function(x00){evaluateLogConDens(x00,x=x,phi,Fhat, IsKnot)[3]}
    ##apply(matrix(x0),1,myf)
    evaluateLogConDens(x0,res=res1,which=3)[,4]
  }
  E.f <- intFfn(x,phi,Fhat)
  {##get phiP=phi prime= deriv of phi; left and right derivs are L and R.
    phiK <- phi[KK]
    slopes <- diff(phiK)/diff(x[KK])
    phiPR.f <- stepfun(x=x[KK], c(slopes[1],slopes,tail(slopes,1)), right=FALSE,f=0)
    phiPL.f <- stepfun(x=x[KK], c(slopes[1],slopes,tail(slopes,1)), right=TRUE,f=1)
    phiPL <- as.vector(phiPL.f(x))
    phiPR <- as.vector(phiPR.f(x))
  }
  ## if (!is.null(logfile) && !is.na(logfile)) sink(NULL) ## just use sink() outsie
  if (print){
    print("ASLC: returning")
  }
  return(c(res1,
           list(phi.f=phi.f, fhat.f=fhat.f, Fhat.f=Fhat.f, E.f=E.f,
           phiPL=phiPL,phiPR=phiPR,
           phiPL.f=phiPL.f,phiPR.f=phiPR.f)))
}




## Haven't written a version of this function for mode (e.g. a
## "logConDens.mode" function) because would have to appropriately change x's
## to z's, etc.

logConDens <- function (x, xgrid = NULL, smoothed = TRUE, print = FALSE, gam = NULL, 
                        xs = NULL, prec=10^-10) 
{
  res1 <- activeSetLogCon(x, xgrid = xgrid, print = print, prec=prec)
  if (identical(smoothed, FALSE)) {
    res <- c(res1)
  }
  if (identical(smoothed, TRUE)) {
    if (identical(xs, NULL)) {
      r <- diff(range(x))
      xs <- seq(min(x) - 0.1 * r, max(x) + 0.1 * r, length = 500)
    }
    smo <- evaluateLogConDens(xs, res1, which = 4:5, gam = gam, 
                              print = print)
    f.smoothed <- smo[, "smooth.density"]
    F.smoothed <- smo[, "smooth.CDF"]
    mode <- xs[f.smoothed == max(f.smoothed)]
    res2 <- list(f.smoothed = f.smoothed, F.smoothed = F.smoothed, 
                 gam = gam, xs = xs, mode = mode)
    res <- c(res1, res2)
  }
  res$smoothed <- smoothed
  class(res) <- "dlc"
  return(res)
}
