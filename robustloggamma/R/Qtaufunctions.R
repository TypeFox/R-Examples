#############################################################
#	All functions in this file are copyrighted
#	Author: C. Agostinelli, A. Marazzi,
#               V.J. Yohai and A. Randriamiharisoa
#	Maintainer e-mail: claudio@unive.it
#	Date: January, 1, 2013
#	Version: 0.1
#	Copyright (C) 2013 C. Agostinelli A. Marazzi,
#                  V.J. Yohai and A. Randriamiharisoa
#############################################################

# Biweight functions

rhoBW <- function(x,k){
k2  <- k*k; k4 <- k2*k2; k6 <- k2*k4
x2  <- x*x; x4 <- x2*x2; x6 <- x4*x2
(3*x2/k2-3*x4/k4+x6/k6)*(abs(x)<k)+(abs(x)>=k)}

psiBW <- function(x,k){
(6/k)*(x/k)*(1-(x/k)^2)^2*(abs(x)<k)}

pspBW <- function(x,k){
k2  <- k*k; k4 <- k2*k2; k6 <- k2*k4
x2  <- x*x; x4 <- x2*x2
(6/k2-36*x2/k4+30*x4/k6)*(abs(x)<k)}

rhoBWdnorm <- function(x,k) {rhoBW(x,k)*dnorm(x)}

#------------------------------------------------------------------------------

MscaleW <- function(u,w,b1,c1,tol){
h  <- 1; it <- 0
s0 <- median(abs(u*w))/.6745
if (s0 > tol) {
 while((h > tol) & (it < 50)){
  it <- it+1
  s1 <- (s0^2)*mean(rhoBW((u*w/s0),c1)) / b1
  s1 <- s1^(1/2)
  h  <- abs(s1-s0)/s0
  s0 <- s1} }
s0}

TauscaleW <- function(u,w,b1,c1,b2,c2,tol) {
  tau <- tol
  s0  <- MscaleW(u,w,b1,c1,tol)
  if (s0 > tol)
    tau <- sqrt(s0^2*mean(rhoBW(u*w/s0,c2)) / b2 )
####  cat(tau, '\n')
  return(tau)
}

#------------------------------------------------------------------------------

RegtauW <- function(x,y,w,b1,c1,b2,c2,N,tol=1e-11,seed=567){
# RegtauW.f is the corresponding fortran function
n <- length(x)
a <- b <- ta <- 1:N
# set.seed(seed)
for (i in 1:N) {
  tt   <- sample(1:n,2)
  xx   <- x[tt]
  yy   <- y[tt]
  b[i] <- (yy[2]-yy[1])/(xx[2]-xx[1])
  a[i] <- yy[1]-b[i]*xx[1]
  rs   <- y - b[i]*x -a[i]
  cc <- order(abs(rs))
  n1 <- floor(n/2)
  cc <- cc[1:n1]
  xxx <- x[cc]
  yyy <- y[cc]
  mx  <- mean(xxx)
  my  <- mean(yyy)
  b[i] <- sum((xxx-mx)*(yyy-my))/sum((xxx-mx)^2)
  a[i] <- my-b[i]*mx
  rs    <- y - b[i]*x - a[i]
  ta[i] <- TauscaleW(rs,w,b1,c1,b2,c2,tol)
}
io <- (1:N)[ta == min(ta)]
io <- min(io)
ao <- a[io]
bo <- b[io]
to <- ta[io]
list(ao=ao,bo=bo,to=to)}

WQtau <- function(ri,w,lgrid,c1=1.547647,c2=6.08,N,maxit,tolr){
# Weighted Qtau estimate (grid optimization)
# set tolr=100 to skip the refinement
  n  <- length(ri); X <- matrix(0,ncol=2,nrow=n); nl <- length(lgrid)
  sL1 <- sL2 <- aL1 <- aL2 <- bL1 <- bL2 <- lgrid
  qq <- list(beta = 0, Tscale = 0, nit=0)
  b1 <- integrate(rhoBWdnorm,lower=-10,upper=10,k=c1)$value
  b2 <- integrate(rhoBWdnorm,lower=-10,upper=10,k=c2)$value
  for (i in 1:nl) {
    lam <- lgrid[i]
    ql  <- qloggamma(p=ppoints(n),lambda=lam)
    z   <- RegtauW.f(ql,ri,w,b1,c1,b2,c2,N)
### z   <- RegtauW(ql,ri,w,b1,c1,b2,c2,N)
    X   <- cbind(1,ql)
    B0  <- c(z$ao,z$bo); s0 <- z$to
    if (s0 > tolr) {
      qq  <- IRLStauW(X,ri,w,inib=B0,iniscale=s0,b1,c1,b2,c2,maxit,tolr)
    } else {
      qq$Tscale <- s0; qq$beta <- B0
    } 
###    cat(i,lgrid[i],z$to,qq$Tscale,qq$nit,"\n")
    sL1[i] <- z$to; sL2[i] <- qq$Tscale
    aL1[i] <- z$ao; aL2[i] <- qq$beta[1]
    bL1[i] <- z$bo; bL2[i] <- qq$beta[2]
  }  
  io1 <- (1:nl)[sL1 == min(sL1)]
##  cat(length(io1), '\n')
  io1 <- min(io1)
  lam.est1 <- lgrid[io1]
  s.est1   <- sL1[io1]
  a.est1   <- aL1[io1]
  b.est1   <- bL1[io1]
  io2 <- (1:nl)[sL2 == min(sL2)]
  io2 <- min(io2)
  lam.est2 <- lgrid[io2]
  s.est2   <- sL2[io2]
  a.est2   <- aL2[io2]
  b.est2   <- bL2[io2]
  list(lam1=lam.est1,Stau1=s.est1,mu1=a.est1,sig1=b.est1,
     lam2=lam.est2,Stau2=s.est2,mu2=a.est2,sig2=b.est2,
     sL1=sL1,sL2=sL2)
}

#------------------------------------------------------------------------------

# Refinement functions

IRLStauW <- function(X,y,w,inib,iniscale,b1,c1,b2,c2,maxit,tol){
# Alfio refinement for weightes Qtau
n <- nrow(X); p <- ncol(X); ps1 <- ps2 <- vector("numeric",n)
res  <- y - X %*% inib
res  <- res*w 
if (iniscale == 0) scale <- median(abs(res))/.6745  else  scale <- iniscale
oldbeta <- inib; betadiff <- 2*tol; iter <- 0
while ((betadiff > tol) && (iter < maxit)) {
  scale     <- sqrt( scale^2 * mean(rhoBW(res/scale,c1)) / b1 )
  scaledres <- res/scale
  Wn.num    <- sum(2*rhoBW(scaledres,c2)-psiBW(scaledres,c2)*scaledres)
  Wn.den    <- sum(psiBW(scaledres,c1)*scaledres)
  Wn        <- Wn.num / Wn.den
  zero      <- scaledres == 0; notz <- scaledres != 0    
  ps1[notz] <- psiBW(scaledres[notz],c1) / scaledres[notz]
  ps1[zero] <- pspBW(scaledres[zero],c1)
  ps2[notz] <- psiBW(scaledres[notz],c2) / scaledres[notz]
  ps2[zero] <- pspBW(scaledres[zero],c2)
  weights   <- (Wn * ps1 + ps2)*w*w
  sqweights <- sqrt(weights)
  xw        <- X * as.vector(sqweights)
  yw        <- y * sqweights
  newbeta   <- qr.coef(qr(xw),yw)
  if (any(!is.finite(newbeta))) {newbeta <- inib; scale <- iniscale; break}
  betadiff  <- sqrt(sum((oldbeta - newbeta)^2))
  res       <- y - X %*% newbeta
  res       <- res*w
  oldbeta   <- newbeta
  iter      <- iter + 1}
if (scale==0) tau <- 0 else tau <- sqrt(scale^2*mean(rhoBW(res/scale,c2)) / b2 )
list(beta=newbeta,Mscale=scale,Tscale=tau,nit=iter)}

#------------------------------------------------------------------------------
## Spostata in Eta.R
## Exp.response <- function(lambda,mu,sigma){
## # expected response in original scale
## zero <- 0.0001
## nm   <- 100000
## if (lambda > zero) {
##   alpha <- 1/lambda^2
##   ca    <- lambda/sigma
##   lMu   <- mu - log(alpha)/ca + lgamma(alpha+1/ca)-lgamma(alpha)
##   Mu    <- exp(lMu) }
## if (abs(lambda) <= zero) Mu <- exp(mu + sigma^2/2)
## if (lambda < -zero) {
##   ql  <- qloggamma(p=ppoints(nm),lambda=lambda)
##   Mu  <- mean(exp(mu+sigma*ql)) }
## Mu}

WQtauboot <- function(ri,w,lambda,step=0.5,lgrid,c1=1.547647,c2=6.08,N,maxit,tolr) {
# Weighted Qtau estimate (for bootstrap)
# set tolr=100 to skip the refinement
  n  <- length(ri); X <- matrix(0,ncol=2,nrow=n)
  b1 <- integrate(rhoBWdnorm,lower=-10,upper=10,k=c1)$value
  b2 <- integrate(rhoBWdnorm,lower=-10,upper=10,k=c2)$value
  temp <- function(lambda) {
    ql  <- qloggamma(p=ppoints(n),lambda=lambda)
##    cat (lambda, 'in temp\n')
    z   <- RegtauW.f(ql,ri,w,b1,c1,b2,c2,N)
    X   <- cbind(1,ql)
    B0  <- c(z$ao,z$bo); s0 <- z$to
    qq  <- IRLStauW(X,ri,w,inib=B0,iniscale=s0,b1,c1,b2,c2,maxit,tolr)
    return(qq)
  }
  qqlambda <- temp(lambda)
  qqlambdadx <- temp(lambda+step)
  qqlambdasx <- temp(lambda-step)
  lambdadx <- lambdasx <- lambda
  tscale <- c(qqlambdasx$Tscale,qqlambda$Tscale,qqlambdadx$Tscale)
## DX
  k <- 1  
  while (qqlambda$Tscale > qqlambdadx$Tscale & lambdadx < max(lgrid)) {
    k <- k + 1
    lambdadx <- lambda+k*step
    qqlambdadx <- temp(lambdadx)
    tscale <- c(tscale, qqlambdadx$Tscale)
  }
  kdx <- k
## SX
  k <- 1  
  while (qqlambda$Tscale > qqlambdasx$Tscale & lambdasx > min(lgrid)) {
    k <- k + 1
    lambdasx <- lambda-k*step
    qqlambdasx <- temp(lambdasx)
    tscale <- c(qqlambdasx$Tscale, tscale)    
  }
  ksx <- k
  lambda <- c(lambda-(ksx:1)*step,lambda,lambda+(1:kdx)*step)
  nl <- length(lambda)
  pos <- which.min(tscale)
  lambdasx <- lambda[max(1,pos-1)]
  lambdadx <- lambda[min(nl,pos+1)]
  
##  cat(lambdasx, lambdadx, pos, nl,"\n")

##  lambdasx1 <- max(lambdasx, min(lgrid))
##  lambdadx1 <- min(lambdadx, max(lgrid))

##  lambdasx <- min(lambdasx1, lambdadx1)
##  lambdadx <- max(lambdasx1, lambdadx1)

##  cat(lambdasx, lambdadx, pos, nl,"\n")

  lgrid <- seq(lambdasx, lambdadx, by=0.05)
  nl <- length(lgrid)
  sL1 <- sL2 <- aL1 <- aL2 <- bL1 <- bL2 <- lgrid
  qq <- list(beta = 0, Tscale = 0, nit=0)
  for (i in 1:nl) {
    lam <- lgrid[i]
##    cat (lam, "\n")
    ql <- qloggamma(p=ppoints(n),lambda=lam)
    z <- RegtauW.f(ql,ri,w,b1,c1,b2,c2,N)
    X   <- cbind(1,ql)
    B0  <- c(z$ao,z$bo); s0 <- z$to
    if (s0 > tolr)
      qq  <- IRLStauW(X,ri,w,inib=B0,iniscale=s0,b1,c1,b2,c2,maxit,tolr)
    else {
      qq$Tscale <- s0
      qq$beta <- B0
    } 
### cat(i,lgrid[i],z$to,qq$Tscale,qq$nit,"\n")
    sL1[i] <- z$to; sL2[i] <- qq$Tscale
    aL1[i] <- z$ao; aL2[i] <- qq$beta[1]
    bL1[i] <- z$bo; bL2[i] <- qq$beta[2]
  }  
  io1 <- (1:nl)[sL1 == min(sL1)]
  io1 <- min(io1)
  lam.est1 <- lgrid[io1]
  s.est1   <- sL1[io1]
  a.est1   <- aL1[io1]
  b.est1   <- bL1[io1]
  io2 <- (1:nl)[sL2 == min(sL2)]
  io2 <- min(io2)
  lam.est2 <- lgrid[io2]
  s.est2   <- sL2[io2]
  a.est2   <- aL2[io2]
  b.est2   <- bL2[io2]
  res <- list(lam1=lam.est1,Stau1=s.est1,mu1=a.est1,sig1=b.est1,
    lam2=lam.est2,Stau2=s.est2,mu2=a.est2,sig2=b.est2,sL1=sL1,sL2=sL2)
  return(res)
}

