#############################################################
#	All functions in this file are copyrighted
#	Author: C. Agostinelli, A. Marazzi,
#               V.J. Yohai and A. Randriamiharisoa
#	Maintainer e-mail: claudio@unive.it
#	Date: February, 04, 2015
#	Version: 0.2
#	Copyright (C) 2015 C. Agostinelli A. Marazzi,
#                  V.J. Yohai and A. Randriamiharisoa
#############################################################
# Disparity.WML_functions.s

weights.loggamma <- function(yi.sorted,mu,sigma,lambda,bw=0.3,raf="SCHI2",tau=0.5,nmod=1000){
  data  <- (yi.sorted-mu)/sigma
  Sdat  <- density(data, kernel="gaussian", bw=bw, cut=3,n=512)
  sdat  <- approxfun(x=Sdat$x, y=Sdat$y, rule=2)
  pe    <- sdat(data)
  modl  <- qloggamma(p=ppoints(nmod),lambda=lambda)
  Smod  <- density(modl, kernel="gaussian", bw=bw, cut=3,n=512)
  smod  <- approxfun(x=Smod$x, y=Smod$y, rule=2)
  pm    <- smod(data)
  delta <- pe/pm-1
  w <- pesi(x=delta, raf=raf, tau=tau)
  w
}

pesi <- function(x, raf, tau=0.1) { # x: residui di Pearson
  gkl <- function(x, tau) {
    if (tau!=0) x <- log(tau*x+1)/tau
    return(x)
  }
  pwd <- function(x, tau) {
    if(tau==Inf)
      x <- log(x+1)
    else
      x <- tau*((x + 1)^(1/tau) - 1)
    return(x)
  }
  x  <- ifelse(x < 1e-10, 0, x)
  ww <- switch(raf,         ## park+basu+2003.pdf
    GKL   = gkl(x, tau),  ## lindsay+1994.pdf                 
    PWD   = pwd(x, tau),
    HD    = 2*(sqrt(x + 1) - 1) ,
    NED   = 2 - (2 + x)*exp(-x) ,
    SCHI2 = 1-(x^2/(x^2 +2)) )
  if (raf!='SCHI2')
    ww <- (ww + 1)/(x + 1)
  ww[ww > 1] <- 1
  ww[ww < 0] <- 0
  ww[is.infinite(x)] <- 0
  return(ww)
}

WMLone.loggamma <- function(yi.sorted,mu0,sig0,lam0,bw=0.3,raf="SCHI2",tau=0.5,nmod=1000,step=1,minw=0.04,nexp=1000) {
  yi   <- yi.sorted; err <- 0; tht <- c(mu0,sig0,lam0); rank <- 3; d <- 100
  wi   <- weights.loggamma(yi,mu0,sig0,lam0,bw=bw,raf=raf,tau=tau,nmod=nmod)
  wi[wi < minw] <- 0  
  yq <- qloggamma(ppoints(nexp), mu=mu0, sigma=sig0, lambda=lam0)
  yq[is.infinite(yq)] <- sign(yq[is.infinite(yq)])*max(abs(yq[is.finite(yq)]))
##  J    <- JacobianB(yi,wi,mu0,sig0,lam0)
  J    <- JacobianB(yq,rep(1, length(yq)),mu0,sig0,lam0) ##Use Expected Fisher Information matrix
  U    <- UscoreB(yi,wi,mu0,sig0,lam0)
  rank <- try(qr(J)$rank)
  if (rank == 3) {
    eJ <- eigen(J)
    if (abs(eJ$values[1]/eJ$values[3]) > d) {
      J <- eJ$vectors%*%diag(eJ$values + (eJ$values[1] - eJ$values[3]*d)/(d - 1))%*%t(eJ$vectors)
    } 
    JI  <- qr.solve(J)
    tht <- tht-step*JI%*%U
  } else
    err = 1
  res <- list(mu=tht[1],sig=tht[2],lam=tht[3],weights=wi,error=err)
  return(res)
}

#reparam.loggamma <- list(gam=function(sigma) sqrt(sigma),
#                         gaminv=function(gam) gam^2,
#                         delta=function(sigma) 2*sqrt(sigma))

WMLone.reparam.loggamma <- function(yi.sorted,mu0,sig0,lam0,bw=0.3,raf="SCHI2",tau=0.5,nmod=1000,step=1,minw=0.04,nexp=1000,reparam) {
  gam0 <- reparam$gam(sig0)
  yi   <- yi.sorted; err <- 0; tht <- c(mu0,gam0,lam0); rank <- 3; d <- 100
  wi   <- weights.loggamma(yi,mu0,sig0,lam0,bw=bw,raf=raf,tau=tau,nmod=nmod)
  wi[wi < minw] <- 0
  yq <- qloggamma(ppoints(nexp), mu=mu0, sigma=sig0, lambda=lam0)
  yq[is.infinite(yq)] <- sign(yq[is.infinite(yq)])*max(abs(yq[is.finite(yq)]))  
##  J    <- JacobianGB(yi,wi,mu0,sig0,lam0,reparam$delta)
  J    <- JacobianGB(yq,rep(1,length(yq)),mu0,sig0,lam0,reparam$delta) ##Use Expected Fisher Information matrix
  U    <- UscoreGB(yi,wi,mu0,sig0,lam0,reparam$delta)
  rank <- try(qr(J)$rank)
  if (rank == 3) {
    eJ <- eigen(J)
    if (abs(eJ$values[1]/eJ$values[3]) > d) {
      J <- eJ$vectors%*%diag(eJ$values + (eJ$values[1] - eJ$values[3]*d)/(d - 1))%*%t(eJ$vectors)
    }     
    JI  <- qr.solve(J)
    tht <- tht-step*JI%*%U
  } else
    err = 1
  res <- list(mu=tht[1],sig=reparam$gaminv(tht[2]),lam=tht[3],weights=wi,error=err)
  return(res)
}

Disparity.WML.loggamma <- function(y,mu0,sig0,lam0,lam.low,lam.sup,tol=0.001,maxit=100,lstep=TRUE,sigstep=TRUE,bw=0.3,raf="SCHI2",tau=0.5,nmod=1000,minw=0.04) {
# solves disparity based WML equations for lambda in (lam.low,lam.sup); mu0, sig0,lam0 are the initial values
  sig.low <- 1e-3*sig0; sig.sup <- 10000; nit <- 1
  dsig <- dlam <- dmu <- 0; p <- 0; n <- length(y)
  sig <- sig0; lam <- lam0; mu <- mu0; dd <- aa <- 0
  repeat {
    wi <- weights.loggamma(y,mu,sig,lam,bw=bw,raf=raf,tau=tau,nmod=nmod)
    wi[wi < minw] <- 0  
# scale step
    sig1 <- sig
    if (sigstep) {
      z <- Nwtsig.loggamma(sig1,mu,lam,y,wi,dd,tol,maxit=maxit)$sig
      if (!is.na(z))
        sig <- z
      else {
        z   <- Rgfsig.loggamma(sig1,mu,lam,y,wi,dd,tol=0.001,maxit=50,sig.low=sig.low,sig.sup=sig.sup)
        sig <- z$sig
        if (sig == -sig1)
          return(list(mu=mu0,sig=sig0,lam=lam0,nit=NA))
      }
    }
    dsig <- sig-sig1
# location step
    mu1     <- mu
    rs      <- (y-mu1)/sig
    di      <- rep(1,n)
    cnd     <- rs!=0
    di[cnd] <- -csiLG(rs[cnd],lam)/rs[cnd]
    mu      <- mean( wi*di*y )/mean( wi*di )
    dmu     <- mu - mu1
# shape step
    lam1 <- lam
    if (lstep) {
      z <- Nwtlam.loggamma(sig,mu,lam1,y,wi,aa,lam.low,lam.sup,tol,maxit=maxit)$lam
      if (!is.na(z))
        lam <- z
      else {
        z <- Rgflam.loggamma(sig,mu,lam1,y,wi,aa,tol=0.001,maxit=50,lam.low,lam.sup)
        lam <- z$lam      
        if (lam == -lam1)
          return(list(mu=mu0,sig=sig0,lam=lam0,nit=NA))
      }
    }
    dlam <- lam-lam1
# cat(nit,lam,mu,sig,"\n")
    if (nit==maxit)
      cat("WML: nit=maxit","\n")
    if (nit==maxit | abs(dmu) < tol  &  abs(dsig) < tol & abs(dlam) < tol ) break
    nit <- nit+1
  }
  list(mu=mu,sig=sig,lam=lam,nit=nit,weights=wi)
}

## if (FALSE) { ### SPOSTATE NEL FILE SCORE&DERIVATIVES.R

## #------------------------------------------------------------------------------

## eta1 <- function(u,lambda){
## # fp.over.f.loggamma(u,lambda)
## csiLG(u,lambda)} 

## eta2 <- function(u,lambda){
## # fp.over.f.loggamma(u,lambda)*u + 1
## csiLG(u,lambda)*u + 1} 

## eta3 <- function(u,lambda){
## # l.score.loggamma(u,lambda)
## psiLG(u,lambda)}

## eta11f <- function(u,lambda){(eta1(u,lambda))^2*dloggamma(x=u, lambda=lambda)}
## eta22f <- function(u,lambda){(eta2(u,lambda))^2*dloggamma(x=u, lambda=lambda)}
## eta33f <- function(u,lambda){(eta3(u,lambda))^2*dloggamma(x=u, lambda=lambda)}
## eta12f <- function(u,lambda){(eta1(u,lambda))*(eta2(u,lambda))*dloggamma(x=u, lambda=lambda)}
## eta13f <- function(u,lambda){(eta1(u,lambda))*(eta3(u,lambda))*dloggamma(x=u, lambda=lambda)}
## eta23f <- function(u,lambda){(eta2(u,lambda))*(eta3(u,lambda))*dloggamma(x=u, lambda=lambda)}

## ###### MESSA NEL FILE ACOV.R E MODIFICATA

## ACOV.ML.loggamma <- function(lambda,lower=-50,upper=50){
## # Asymptotic covariance matrix of ML estimates for standard loggamma: mu=0, sigma=1
## J11 = integrate( eta11f,lower=lower,upper=upper,lambda=lambda)$value
## J22 = integrate( eta22f,lower=lower,upper=upper,lambda=lambda)$value
## J33 = integrate( eta33f,lower=lower,upper=upper,lambda=lambda)$value
## J12 = integrate( eta12f,lower=lower,upper=upper,lambda=lambda)$value
## J13 = integrate( eta13f,lower=lower,upper=upper,lambda=lambda)$value
## J23 = integrate( eta23f,lower=lower,upper=upper,lambda=lambda)$value
## Jf  = matrix(c(J11,J12,J13,J12,J22,J23,J13,J23,J33),nrow=3,byrow=TRUE)
## list(Fisher.Info=Jf, ACOV=solve(Jf))}

## }
