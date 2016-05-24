#############################################################
# 
#	Internal functions and wrappers from Fortran
#	Author: Claudio Agostinelli and Victor J. Yohai
#	E-mail: claudio@unive.it
#	Date: June, 22, 2014
#	Version: 0.1
#
#	Copyright (C) 2014 Claudio Agostinelli
#                      and Victor J. Yohai
#
#############################################################

vcrobresid <- function(x, y, beta) {
  res <- y - xprod(x,beta)
  return(res)
}

rssr <- function(resid, Sigma) {
  p <- nrow(resid)
  n <- ncol(resid)
  JL <- p*(p-1)/2

  rssr <- .Fortran("rsstarr",
    as.matrix(resid),
    as.integer(p), 
    as.integer(n),
    as.integer(JL),
    as.matrix(Sigma),
    rssr=matrix(0, nrow=JL, ncol=n),
    package="varComprob")$rssr
  return(rssr)
}

rsr <- function(resid, Sigma) {
  p <- nrow(resid)
  n <- ncol(resid)
  JL <- p*(p-1)/2

  rsr <- .Fortran("rsr",
    as.matrix(resid),
    as.integer(p), 
    as.integer(n),
    as.integer(JL),
    as.matrix(Sigma),
    rsr=matrix(0, nrow=JL, ncol=n),
    package="varComprob")$rsr
  return(rsr)
}

xssx <- function(x, Sigma) {
  dx <- dim(x)
  p <- dx[1]
  n <- dx[2]
  k <- dx[3]
  JL <- p*(p-1)/2
  
  xssx <- .Fortran("xsstarx",
    as.array(x),
    as.integer(p), 
    as.integer(n),
    as.integer(k),                   
    as.integer(JL),
    as.matrix(Sigma),
    xssx=array(0.0, dim=c(JL,n,k,k)),
    package="varComprob")$xssx
  return(xssx)
}

xssy <- function(x, y, Sigma) {
  dx <- dim(x)
  p <- dx[1]
  n <- dx[2]
  k <- dx[3]
  JL <- p*(p-1)/2
  
  xssy <- .Fortran("xsstary",
    as.array(x),
    as.matrix(y),
    as.integer(p), 
    as.integer(n),
    as.integer(k),                   
    as.integer(JL),
    as.matrix(Sigma),
    xssy=array(0.0, dim=c(JL,n,k)),
    package="varComprob")$xssy
  return(xssy)
}

rvr <- function(resid, V) {
  p <- nrow(resid)
  n <- ncol(resid)
  JL <- p*(p-1)/2
  R <- dim(V)[3]

  rvr <- .Fortran("rvr",
    as.matrix(resid),
    as.integer(p), 
    as.integer(n),
    as.integer(JL),
    as.array(V),
    as.integer(R),
    rvr=array(0.0, dim=c(JL,n,R)),
    package="varComprob")$rvr
  return(rvr)
}

sdet <- function(Sigma) {
  p <- nrow(Sigma)
  JL <- p*(p-1)/2

  sigmadet <- .Fortran("sdet",
    as.matrix(Sigma),
    as.integer(p),
    as.integer(JL),
    sigmadet=double(JL),
    package="varComprob")$sigmadet
  return(sigmadet)
}

## GOAL FUNCTIONS
GoalCompositeS <- function(x, Y, X, scales, V, controllo) {
  ## x: c(beta,gamma)
  ## Y: matrix. dim(y)=c(p,n) 
  ## X: array. dim(x)=c(p,n,k)
  ## V: array. dim(V)=c(p,p,R)
  dX <- dim(X)
  p <- dX[1]
  n <- dX[2]
  k <- dX[3]  
  dV <- dim(V)
  R <- dV[3]
  beta <- x[1:k]
  gamma <- x[(k+1):(k+R)]
  resid <- as.matrix(vcrobresid(y=Y, x=X, beta=beta))  
  Sigma <- Vprod(V=V, gamma=gamma) ## V0 is added automatically in Vprod
  if (any(eigen(Sigma)$values <= 0))
    return(2*sum(scales)+1)
  RR <- as.matrix(rssr(resid=resid, Sigma=Sigma))
  scales <- doSsteppw(RR=RR, scale=scales, bb=controllo$bb, cc=controllo$tuning.chi, psi=controllo$psi, tol=controllo$rel.tol.scale, verbose=(controllo$trace.lev>2))
  S <- sum(scales)
  return(S)
}

GoalClassicS <- function(x, Y, X, scale, V, controllo, nsize=NULL) {
  ## x: c(beta,gamma)
  ## Y: matrix. dim(y)=c(p,n) 
  ## X: array. dim(x)=c(p,n,k)
  ## V: array. dim(V)=c(p,p,R)
  dX <- dim(X)
  p <- dX[1]
  n <- dX[2]
  k <- dX[3]  
  dV <- dim(V)
  R <- dV[3]
  beta <- x[1:k]
  gamma <- x[(k+1):(k+R)]
  resid <- as.matrix(vcrobresid(y=Y, x=X, beta=beta))  
  Sigma <- Vprod(V=V, gamma=gamma) ## V0 is added automatically
  if (any(eigen(Sigma)$values <= 0))
    return(2*scale+1)
  Sigmastar <- Sigma/det(Sigma)^(1/nrow(Sigma))
  Sigmastarinv <- solve(Sigmastar)
  RR <- rep(0, n)
  for (i in 1:n)
    RR[i] <- drop(resid[,i]%*%Sigmastarinv%*%resid[,i])
  if (controllo$psi!="rocke")
    scale <- doSstep(m=RR, scale=scale, bb=controllo$bb, cc=controllo$tuning.chi, psi=controllo$psi, tol=controllo$rel.tol.scale, verbose=(controllo$trace.lev>2))
  else {
    if (!is.null(nsize))
      RR <- rep(RR,nsize)
    scale <- doSsteprocke(m=RR, scale=scale, bb=controllo$bb, p=p, arp=controllo$arp.chi, tol=controllo$rel.tol.scale/1000, verbose=(controllo$trace.lev>2))
  }
  return(scale)
}

GoalCompositeMM <- function(x, Y, X, scales, V, controllo) {
  ## only one observation at time!
  ## x: c(beta,gamma)
  ## Y: matrix. dim(y)=c(p,1) 
  ## X: array. dim(x)=c(p,1,k)
  ## V: array. dim(V)=c(p,p,R)
  dX <- dim(X)
  p <- dX[1]
  JL <- p*(p-1)/2
  n <- dX[2]
  k <- dX[3]  
  dV <- dim(V)
  R <- dV[3]
  beta <- x[1:k]
  gamma <- x[(k+1):(k+R)]
  resid <- as.matrix(vcrobresid(y=Y, x=X, beta=beta))  
  Sigma <- Vprod(V=V, gamma=gamma) ## V0 is added automatically in Vprod
  if (any(eigen(Sigma)$values <= 0))
    return(2*sum(scales)+1)
  RR <- as.matrix(rssr(resid=resid, Sigma=Sigma))
##  S <- sum(scales*Mchi(sqrt(RR/scales), cc=controllo$tuning.psi, psi=controllo$psi))
  S <- sum(sweep(x=matrix(Mchi(sqrt(sweep(x=RR,MARGIN=1,STATS=scales,FUN="/")), cc=controllo$tuning.psi, psi=controllo$psi), nrow=JL, ncol=n, byrow=FALSE), MARGIN=1,STATS=scales,FUN="*"))  
  return(S)
}

GoalClassicMM <- function(x, Y, X, scale, V, controllo) {
  ## x: c(beta,gamma)
  ## Y: matrix. dim(y)=c(p,n) 
  ## X: array. dim(x)=c(p,n,k)
  ## V: array. dim(V)=c(p,p,R)
  dX <- dim(X)
  p <- dX[1]
  n <- dX[2]
  k <- dX[3]  
  dV <- dim(V)
  R <- dV[3]
  beta <- x[1:k]
  gamma <- x[(k+1):(k+R)]
  resid <- as.matrix(vcrobresid(y=Y, x=X, beta=beta))  
  Sigma <- Vprod(V=V, gamma=gamma) ## V0 is added automatically
  if (any(eigen(Sigma)$values <= 0))
    return(2*scale+1)
  Sigmastar <- Sigma/det(Sigma)^(1/nrow(Sigma))
  Sigmastarinv <- solve(Sigmastar)
  RR <- rep(0, n)
  for (i in 1:n)
    RR[i] <- drop(resid[,i]%*%Sigmastarinv%*%resid[,i])
  if (controllo$psi!="rocke")
    S <- sum(Mchi(x=sqrt(RR/scale), cc=controllo$tuning.psi, psi=controllo$psi))
  else
    S <- sum(rho.rk2.f(x=(RR/scale), p=p, alpha=controllo$arp.psi))
  return(S)
}

GoalCompositeTau <- function(x, Y, X, scales, V, controllo) {
  ## x: c(beta,gamma)
  ## Y: matrix. dim(y)=c(p,n) 
  ## X: array. dim(x)=c(p,n,k)
  ## V: array. dim(V)=c(p,p,R)
  dX <- dim(X)
  p <- dX[1]
  n <- dX[2]
  k <- dX[3]  
  dV <- dim(V)
  R <- dV[3]
  beta <- x[1:k]
  gamma <- x[(k+1):(k+R)]
  resid <- as.matrix(vcrobresid(y=Y, x=X, beta=beta))  
  Sigma <- Vprod(V=V, gamma=gamma) ## V0 is added automatically in Vprod
  if (any(eigen(Sigma)$values <= 0))
    return(2*sum(scales)+1)
  RR <- as.matrix(rssr(resid=resid, Sigma=Sigma))
  scales <- doSsteppw(RR=RR, scale=scales, bb=controllo$bb, cc=controllo$tuning.chi, psi=controllo$psi, tol=controllo$rel.tol.scale, verbose=(controllo$trace.lev>2))
  scales <- doTausteppw(RR=RR, scale=scales, cc=controllo$tuning.psi, psi=controllo$psi)    
  T <- sum(scales)
  return(T)
}

GoalClassicTau <- function(x, Y, X, scale, V, controllo) {
  .NotYetImplemented()
}

## BETAS
doBetastep <- function(Wdot, XX, XY) {
  dx <- dim(XX)
  JL <- dx[1]
  N <- dx[2]
  k <- dx[3]
  num <- matrix(0, nrow=k, ncol=1)
  den <- matrix(0, nrow=k, ncol=k)
  for (n in 1:N) {
    for (jl in 1:JL) {
      num <- num + Wdot[jl,n]*XY[jl,n,]
      den <- den + Wdot[jl,n]*XX[jl,n,,]
    }
  }
  beta <- solve(den)%*%num
  return(beta)
}

doGammaCompositeSGoal <- function(x, resid, scales, V, controllo) {
  ## x: gamma
  dV <- dim(V)
  p <- dV[1]
  JL <- p*(p-1)/2
  n <- ncol(resid)
  Sigma <- Vprod(V=V, gamma=x) ## V0 is added automatically
  if (any(eigen(Sigma)$values <= 0))
    return(2*sum(scales)+1)
  RR <- rssr(resid=resid, Sigma=Sigma)
  scales <- doSsteppw(RR=RR, scale=scales, bb=controllo$bb, cc=controllo$tuning.chi, psi=controllo$psi, tol=controllo$rel.tol.scale, verbose=(controllo$trace.lev>2))
  S <- sum(scales)
  return(S)
}

## GAMMAS Algoritmo Optim
doGammaCompositeSstep <- function(gamma, resid, scales, V, control, ...) {
  lower <- rep(control$lower, length.out=length(gamma))
  upper <- rep(control$upper, length.out=length(gamma))
  res <- optim(par=gamma, fn=doGammaCompositeSGoal, method="L-BFGS-B", lower=lower, upper=upper, resid=resid, scales=scales, V=V, controllo=control, ...)
  return(res$par)
}

doGammaCompositeMMGoal <- function(x, resid, scales, V, Mmax, controllo) {
  ## x: gamma
  dV <- dim(V)
  p <- dV[1]
  JL <- p*(p-1)/2
  n <- ncol(resid)
  Sigma <- Vprod(V=V, gamma=x) ## V0 is added automatically
  if (any(eigen(Sigma)$values <= 0))
    return(Mmax)
  RR <- rssr(resid=resid, Sigma=Sigma)
  S <- sum(sweep(x=matrix(Mchi(sqrt(sweep(x=RR,MARGIN=1,STATS=scales,FUN="/")), cc=controllo$tuning.psi, psi=controllo$psi), nrow=JL, ncol=n, byrow=FALSE), MARGIN=1,STATS=scales,FUN="*"))
  return(S)
}

doGammaCompositeMMstep <- function(gamma, resid, scales, V, Mmax, control, ...) {
  lower <- rep(control$lower, length.out=length(gamma))
  upper <- rep(control$upper, length.out=length(gamma))
  res <- optim(par=gamma, fn=doGammaCompositeMMGoal, method="L-BFGS-B", lower=lower, upper=upper, resid=resid, scales=scales, V=V, Mmax=Mmax, controllo=control, ...)
  return(res$par)
}

doGammaCompositeTauGoal <- function(x, resid, scales, V, Tmax, controllo) {
  ## x: gamma
  dV <- dim(V)
  p <- dV[1]
  JL <- p*(p-1)/2
  n <- ncol(resid)
  Sigma <- Vprod(V=V, gamma=x) ## V0 is added automatically
  if (any(eigen(Sigma)$values <= 0))
    return(Tmax)
  RR <- rssr(resid=resid, Sigma=Sigma)
###    browser(expr=any(is.nan(RR)))
  scales <- doSsteppw(RR=RR, scale=scales, bb=controllo$bb, cc=controllo$tuning.chi, psi=controllo$psi, tol=controllo$rel.tol.scale, verbose=(controllo$trace.lev>2))
  if (controllo$psi=="optimal")
    scales <- scales*controllo$tuning.chi^2
  T <- doTausteppw(RR=RR, scale=scales, cc=controllo$tuning.psi, psi=controllo$psi)
  return(T)
}

## GAMMAS Algoritmo for Tau Optim
doGammaCompositeTaustep <- function(gamma, resid, scales, V, Tmax, control, ...) {
  lower <- rep(control$lower, length.out=length(gamma))
  upper <- rep(control$upper, length.out=length(gamma))    
  res <- optim(par=gamma, fn=doGammaCompositeTauGoal, method="L-BFGS-B", lower=lower, upper=upper, resid=resid, scales=scales, V=V, Tmax=Tmax, controllo=control, ...)
  if (control$trace.lev > 1) {
    cat('Gamma Step\n')
    cat('gamma: ', res$par, '\n')
    cat('value of the function: ', res$value, '\n')
  }
  return(res$par)
}

#### GAMMA for CLASSIC VERSIONS
doGammaClassicSGoal <- function(x, resid, scale, V, controllo) {
    ## x: gamma
    dV <- dim(V)
    p <- dV[1]
    n <- ncol(resid)
    Sigma <- Vprod(V=V, gamma=x) ## V0 is added automatically
    if (any(eigen(Sigma)$values <= 0))
      return(2*scale+1)
    Sigmastar <- Sigma/det(Sigma)^(1/nrow(Sigma))
    Sigmastarinv <- solve(Sigmastar)
    RR <- rep(0, n)
    for (i in 1:n)
      RR[i] <- drop(resid[,i]%*%Sigmastarinv%*%resid[,i])
    if (controllo$psi!="rocke")  
      scale <- doSstep(m=RR, scale=10, bb=controllo$bb, cc=controllo$tuning.chi, psi=controllo$psi, tol=controllo$rel.tol.scale, verbose=(controllo$trace.lev>2))
    else
      scale <- doSsteprocke(m=RR, scale=10, bb=controllo$bb, p=p, arp=controllo$arp.chi, tol=controllo$rel.tol.scale/1000, verbose=(controllo$trace.lev>2))
    return(scale)
}

doGammaClassicSstep <- function(gamma, resid, scale, V, control, ...) {
  lower <- rep(control$lower, length.out=length(gamma))
  upper <- rep(control$upper, length.out=length(gamma))
  res <- optim(par=gamma, fn=doGammaClassicSGoal, method="L-BFGS-B", lower=lower, upper=upper, resid=resid, scale=scale, V=V, controllo=control, ...)
  return(res$par)
}

doGammaClassicMMGoal <- function(x, resid, scale, V, Mmax, controllo) {
  ## x: gamma
  dV <- dim(V)
  p <- dV[1]
  n <- ncol(resid)
  Sigma <- Vprod(V=V, gamma=x) ## V0 is added automatically
  if (any(eigen(Sigma)$values <= 0))
    return(Mmax)
  Sigmastar <- Sigma/det(Sigma)^(1/nrow(Sigma))
  Sigmastarinv <- solve(Sigmastar)
  RR <- rep(0, n)
  for (i in 1:n)
    RR[i] <- drop(resid[,i]%*%Sigmastarinv%*%resid[,i])
  if (controllo$psi!="rocke")
    S <- sum(Mchi(x=sqrt(RR/scale), cc=controllo$tuning.psi, psi=controllo$psi))
  else
    S <- sum(rho.rk2.f(x=(RR/scale), p=p, alpha=controllo$arp.psi))
  return(S)
}

doGammaClassicMMstep <- function(gamma, resid, scale, V, Mmax, control, ...) {
  lower <- rep(control$lower, length.out=length(gamma))
  upper <- rep(control$upper, length.out=length(gamma))
  res <- optim(par=gamma, fn=doGammaClassicMMGoal, method="L-BFGS-B", lower=lower, upper=upper, resid=resid, scale=scale, V=V, Mmax=Mmax, controllo=control, ...)
  return(res$par)
}

doGammaClassicTauGoal <- function(x, resid, scale, V, Tmax, controllo) {
  ## x: gamma
  dV <- dim(V)
  p <- dV[1]
  n <- ncol(resid)
  Sigma <- Vprod(V=V, gamma=x) ## V0 is added automatically
  if (any(eigen(Sigma)$values <= 0))
    return(Tmax)  
  Sigmastar <- Sigma/det(Sigma)^(1/nrow(Sigma))
  Sigmastarinv <- solve(Sigmastar)
  RR <- rep(0, n)
  for (i in 1:n)
    RR[i] <- drop(resid[,i]%*%Sigmastarinv%*%resid[,i])
  if (controllo$psi!="rocke") {
    scale <- doSstep(m=RR, scale=10, bb=controllo$bb, cc=controllo$tuning.chi, psi=controllo$psi, tol=controllo$rel.tol.scale, verbose=(controllo$trace.lev>2))
    T1 <- mean(rhostar(m=RR, scale=scale, cc=controllo$tuning.psi,  psi=controllo$psi))
  } else {
    scale <- doSsteprocke(m=RR, scale=10, bb=controllo$bb, p=p, arp=controllo$arp.chi, tol=controllo$rel.tol.scale, verbose=(controllo$trace.lev>2))
    T1 <- mean(rho.rk2.f(x=(RR/scale), p=p, alpha=controllo$arp.psi))
  }
  T <- scale*T1
  return(T)
}

doGammaClassicTaustep <- function(gamma, resid, scale, V, Tmax, control, ...) {
  lower <- rep(control$lower, length.out=length(gamma))
  upper <- rep(control$upper, length.out=length(gamma))
  res <- optim(par=gamma, fn=doGammaClassicTauGoal, method="L-BFGS-B", lower=lower, upper=upper, resid=resid, scale=scale, V=V, Tmax=Tmax, controllo=control, ...)
  if (control$trace.lev > 1) {
    cat('Gamma Step\n')
    cat('gamma: ', res$par, '\n')
    cat('value of the function: ', res$value, '\n')
  }
  return(res$par)
}

##SCALES
doSsteppw <- function(RR, scale, bb, cc, psi, tol=10^(-8), verbose=FALSE) {
  JL <- nrow(RR)
  n <- ncol(RR)
  if (psi=="bisquare")
    psi <- 1L
  else if (psi=="optimal")
    psi <- 3L
  S <- .Fortran("dospw",
      as.matrix(RR),
      as.integer(JL),
      as.integer(n),
      as.double(scale),
      as.double(bb),
      as.double(cc),
      as.integer(psi),                
      as.double(tol),
      package="varComprob")[[4]]
  return(S)
}

doSstep <- function(m, scale, bb, cc, psi, tol=10^(-8), verbose=FALSE) {
  n <- length(m)
  if (psi=="bisquare")
    psi <- 1L
  else if (psi=="optimal")
    psi <- 3L
  scale <- .Fortran("dosstep",
    as.double(m),
    as.integer(n),
    as.double(scale),
    as.double(bb),
    as.double(cc),
    as.integer(psi),
    as.double(tol),
    package="varComprob")[[3]]
  return(scale)
}

doSsteprocke <- function(m, scale, bb, p, arp, tol=10^(-8), verbose=FALSE) {
  n <- length(m)
  dq  <- qchisq(1-arp, p)
  if (p >= n)
    stop("Rocke rho function can not be used when p is greater than n")
  scale <- .Fortran("dosstepr",
    as.double(m),
    as.integer(n),
    as.double(scale),
    as.double(bb),                    
    as.integer(p),
    as.double(dq),
    as.double(tol),
    package="varComprob")[[3]]
  return(scale)
}

###################################
## Functions for the TAU estimate
###################################
rhostar <- function(m, scale, cc, psi) {
  n <- length(m)
  if (psi=="bisquare")
    m <- .Fortran("stukeych",
      as.double(m),
      as.integer(n),
      as.double(scale),
      as.double(cc),
      package="varComprob")[[1]]
  else if (psi=="optimal")
    m <- .Fortran("soptimch",
      as.double(m),
      as.integer(n),
      as.double(scale),
      as.double(cc),
      package="varComprob")[[1]]
##    m <- Mchi(x=sqrt(m/scale), cc=cc, psi=psi)
  return(m)
}

# This is A_{k,jl}
rhostarpw <- function(RR, scale, cc, psi) {
  JL <- nrow(RR)
  n <- ncol(RR)
  if (psi=="bisquare")
    psi <- 1L
  else if (psi=="optimal")
    psi <- 3L
  A <- .Fortran("rhospw",
    as.double(RR),
    as.integer(JL),
    as.integer(n),
    as.double(scale),
    as.double(cc),
    as.integer(psi),
    A=double(JL),
    package="varComprob")$A
  return(A)
}

# This is B_{k,jl}
vcrobweights2pw <- function(RR, scale, cc, psi) {
  JL <- nrow(RR)
  B <- matrix(0, JL, ncol(RR))
  for (jl in 1:JL)
    B[jl,] <- vcrobweights(m=RR[jl,], scale=scale[jl], cc=cc, psi=psi)*RR[jl,]/scale[jl]
  B <- apply(B, 1, mean)
  return(B)
}

doTausteppw <- function(RR, scale, cc, psi) {
  JL <- nrow(RR)
  n <- ncol(RR)
  if (psi=="bisquare")
    psi <- 1L
  else if (psi=="optimal")
    psi <- 3L
  T <- .Fortran("dotstep",
    as.double(RR),
    as.integer(JL),
    as.integer(n),
    as.double(scale),
    as.double(cc),
    as.integer(psi),        
    T=double(1),
    package="varComprob")$T
####    T <- sum(scale*rhostarpw(RR, scale, cc))
  return(T)
}

doTausteppwDet <- function(RR, scale, cc, psi, detS) {
  JL <- nrow(RR)
  n <- ncol(RR)
  if (psi=="bisquare")
    psi <- 1L
  else if (psi=="optimal")
    psi <- 3L
  T <- .Fortran("dotstepd",
    as.double(RR),
    as.integer(JL),
    as.integer(n),
    as.double(scale),
    as.double(cc),
    as.integer(psi),
    as.double(detS),
    T=double(1),
    package="varComprob")$T
####    T <- sum(scale*rhostarpw(RR, scale, cc))
  return(T)
}

##WEIGHTS
vcrobweightsdotpw <- function(W, RR, scale) {
  JL <- nrow(RR)
  N <- ncol(RR)
  Wdot <- matrix(0, nrow=nrow(W), ncol=ncol(W))
  for (jl in 1:JL)
    Wdot[jl,] <- N*W[jl,]*scale[jl]/W[jl,]%*%RR[jl,]
  return(Wdot)
}

vcrobweightspw <- function(RR, scale, cc, psi) {
  JL <- nrow(RR)
  W <- matrix(0, JL, ncol(RR))
  for (jl in 1:JL)
    W[jl,] <- vcrobweights(m=RR[jl,], scale=scale[jl], cc=cc, psi=psi)
  return(W)
}

vcrobweights <- function(m, scale, cc, psi) {
  w <- drop(Mpsi(x=sqrt(m/scale), cc=cc, psi=psi, deriv=0)/sqrt(m/scale))
  w[m==0] <- 1
  return(w)
}

## Rocke weights
vcrobweightsrocke <- function(m, scale, p, arp) {
  w <- w.rk2(x=m/scale, p=p, alpha=arp)  
  w[m==0] <- 1
  return(w)
}

### STANDARD ERRORS FOR COMPOSITE S
VCOV.CompositeS <- function(beta, gamma, scales, y, x, V, control) {
  dX <- dim(x)
  p <- dX[1]
  n <- dX[2]
  k <- dX[3]  
  dV <- dim(V)
  R <- dV[3]
  nr <- k + R
  Var <- matrix(NA, nrow=nr, ncol=n)
  H <- matrix(0, nrow=nr, ncol=nr)
  for (i in 1:n) {
    Var[,i] <- grad(func=GoalCompositeS,  x=c(beta, gamma), Y=y[,i,drop=FALSE], X=x[,i,,drop=FALSE], scales=scales, V=V, controllo=control)
    H <- H + hessian(func=GoalCompositeS,  x=c(beta, gamma), Y=y[,i,drop=FALSE], X=x[,i,,drop=FALSE], scales=scales, V=V, controllo=control)
  }
  Var <- cov(t(Var))
  H <- solve(H)
  vcov <- H%*%Var%*%t(H)*n
}

### STANDARD ERRORS FOR COMPOSITE MM
VCOV.CompositeMM <- function(beta, gamma, scales, y, x, V, control) {
  dX <- dim(x)
  p <- dX[1]
  n <- dX[2]
  k <- dX[3]  
  dV <- dim(V)
  R <- dV[3]
  nr <- k + R
  Var <- matrix(NA, nrow=nr, ncol=n)
  H <- matrix(0, nrow=nr, ncol=nr)
  for (i in 1:n) {
    Var[,i] <- grad(func=GoalCompositeMM,  x=c(beta, gamma), Y=y[,i,drop=FALSE], X=x[,i,,drop=FALSE], scales=scales, V=V, controllo=control)
    H <- H + hessian(func=GoalCompositeMM,  x=c(beta, gamma), Y=y[,i,drop=FALSE], X=x[,i,,drop=FALSE], scales=scales, V=V, controllo=control)
  }
  Var <- cov(t(Var))
  H <- solve(H)
  vcov <- H%*%Var%*%t(H)*n
}


### STANDARD ERRORS FOR COMPOSITE TAU
VCOV.CompositeTau <- function(beta, gamma, scales, y, x, V, control) {
  dX <- dim(x)
  p <- dX[1]
  n <- dX[2]
  k <- dX[3]  
  dV <- dim(V)
  R <- dV[3]
  nr <- k + R
  Var <- matrix(NA, nrow=nr, ncol=n)
  H <- matrix(0, nrow=nr, ncol=nr)
  for (i in 1:n) {
    Var[,i] <- grad(func=GoalCompositeTau,  x=c(beta, gamma), Y=y[,i,drop=FALSE], X=x[,i,,drop=FALSE], scales=scales, V=V, controllo=control)
    H <- H + hessian(func=GoalCompositeTau,  x=c(beta, gamma), Y=y[,i,drop=FALSE], X=x[,i,,drop=FALSE], scales=scales, V=V, controllo=control)
  }
  Var <- cov(t(Var))
  H <- solve(H)
  vcov <- H%*%Var%*%t(H)*n
}

### STANDARD ERRORS FOR CLASSIC S
VCOV.ClassicS <- function(beta, gamma, scale, y, x, V, control) {
  dX <- dim(x)
  p <- dX[1]
  n <- dX[2]
  k <- dX[3]  
  dV <- dim(V)
  R <- dV[3]
  nr <- k + R
  Var <- matrix(NA, nrow=nr, ncol=n)
  H <- matrix(0, nrow=nr, ncol=nr)
  if (control$psi=="rocke")
    nsize <- n
  else
    nsize <- NULL
  for (i in 1:n) {
    Var[,i] <- grad(func=GoalClassicS,  x=c(beta, gamma), Y=y[,i,drop=FALSE], X=x[,i,,drop=FALSE], scale=scale, V=V, controllo=control, nsize=nsize)
    H <- H + hessian(func=GoalClassicS,  x=c(beta, gamma), Y=y[,i,drop=FALSE], X=x[,i,,drop=FALSE], scale=scale, V=V, controllo=control, nsize=nsize)
  }
  Var <- cov(t(Var))
  H <- solve(H)
  vcov <- H%*%Var%*%t(H)*n
}

### STANDARD ERRORS FOR CLASSIC MM
VCOV.ClassicMM <- function(beta, gamma, scale, y, x, V, control) {
  dX <- dim(x)
  p <- dX[1]
  n <- dX[2]
  k <- dX[3]  
  dV <- dim(V)
  R <- dV[3]
  nr <- k + R
  Var <- matrix(NA, nrow=nr, ncol=n)
  H <- matrix(0, nrow=nr, ncol=nr)
  for (i in 1:n) {
    Var[,i] <- grad(func=GoalClassicMM,  x=c(beta, gamma), Y=y[,i,drop=FALSE], X=x[,i,,drop=FALSE], scale=scale, V=V, controllo=control)
    H <- H + hessian(func=GoalClassicMM,  x=c(beta, gamma), Y=y[,i,drop=FALSE], X=x[,i,,drop=FALSE], scale=scale, V=V, controllo=control)
  }
  Var <- cov(t(Var))
  H <- solve(H)
  vcov <- H%*%Var%*%t(H)*n
}


### STANDARD ERRORS FOR CLASSIC TAU
VCOV.ClassicTau <- function(beta, gamma, scale, y, x, V, control) {
  dX <- dim(x)
  p <- dX[1]
  n <- dX[2]
  k <- dX[3]  
  dV <- dim(V)
  R <- dV[3]
  nr <- k + R
  Var <- matrix(NA, nrow=nr, ncol=n)
  H <- matrix(0, nrow=nr, ncol=nr)
  for (i in 1:n) {
    Var[,i] <- grad(func=GoalClassicTau,  x=c(beta, gamma), Y=y[,i,drop=FALSE], X=x[,i,,drop=FALSE], scale=scale, V=V, controllo=control)
    H <- H + hessian(func=GoalClassicTau,  x=c(beta, gamma), Y=y[,i,drop=FALSE], X=x[,i,,drop=FALSE], scale=scale, V=V, controllo=control)
  }
  Var <- cov(t(Var))
  H <- solve(H)
  vcov <- H%*%Var%*%t(H)*n
}

xprod <- function(x, beta) {
  beta <- matrix(beta, ncol=1)
  dx <- dim(x)
  p <- dx[1]
  N <- dx[2]
  k <- dx[3]
  X <- matrix(0, nrow=p, ncol=N)
  for (n in 1:N) {
    X[,n] <- drop(x[,n,]%*%beta)
  }
  return(X)
}

Vprod <- function(V, gamma) {
  dV <- dim(V)
  p <- dV[1]
  R <- dV[3]
  VV <- diag(p) ## V_0
  for (r in 1:R) {
    VV <- VV+V[,,r]*gamma[r]
  }
  return(VV)
}
