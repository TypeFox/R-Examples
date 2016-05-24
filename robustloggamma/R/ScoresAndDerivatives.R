#############################################################
#	All functions in this file are copyrighted
#	Author: C. Agostinelli, A. Marazzi,
#               V.J. Yohai and A. Randriamiharisoa
#	Maintainer e-mail: claudio@unive.it
#	Date: January, 16, 2014
#	Version: 0.1-1
#	Copyright (C) 2014 C. Agostinelli A. Marazzi,
#                  V.J. Yohai and A. Randriamiharisoa
#############################################################

# Extended loggamma distribution
# =======================================

dloggamma <- function(x, mu=0, sigma=1, lambda, log = FALSE) {
# generalized loggamma density
  zero <- 0.0001
  x <- (x-mu)/sigma
  if (abs(lambda) > zero) {
    lam2 <- lambda^(-2)
    res  <- log(abs(lambda))+lam2*log(lam2)+(lam2)*(lambda*x-exp(lambda*x))-lgamma(lam2) - log(sigma)
    if (!log)
      res <- exp(res)
  } else {
    res <- dnorm(x, mean = 0, sd = 1, log = log)/sigma
  }
  return(res)
}

ploggamma <- function(q, mu=0, sigma=1, lambda, lower.tail=TRUE, log.p=FALSE) { 
# generalized loggamma cdf
  zero <- 0.0001
  q <- (q-mu)/sigma
  if (lambda < - zero)
    lower.tail <- !lower.tail
  if (abs(lambda) > zero) {
    alpha <- 1/lambda^2
    res   <- pgamma(alpha*exp(lambda*q), shape=alpha, rate=1, lower.tail=lower.tail, log.p=log.p)
  } else {
    res <- pnorm(q, mean=0, sd=1, lower.tail=lower.tail, log.p=log.p)
  }
  return(res)
}

## ploggamma <- function(q, mu=0, sigma=1, lambda) { 
## # generalized loggamma cdf
##   n <- length(q)
##   if (length(mu)!=1 | length(sigma)!=1 | length(lambda)!=1)
##     stop("Parameters must have length equal to one")
##   res <- .Fortran("ploggamm",
##     as.double(q),
##     as.integer(n),
##     as.double(mu),
##     as.double(sigma),
##     as.double(lambda),
##     p = double(n),
##     PACKAGE="robustloggamma"
##   )
##   return(res$p)
## }

qloggamma <- function(p, mu=0, sigma=1, lambda) { 
# p-quantile of generalized loggamma cdf
  zero <- 0.0001
  if (lambda < -zero) p <- 1-p
  if (abs(lambda) > zero) {
    k   <- 1/lambda^2
    res <- (log(qgamma(p,shape=k,rate=1))+log(lambda^2))/lambda
  } else {
    res <- qnorm(p,mean=0,sd=1)
  }
  res <- res*sigma + mu
  return(res)
}

rloggamma <- function(n, mu=0, sigma=1, lambda) { 
# generates n random generalized loggamma variates
  zero <- 0.0001
  if (abs(lambda) > zero) { 
    alpha <- 1/lambda^2
    ei    <- rgamma(n,shape=alpha,rate=1)
    ui    <- (log(ei)-log(alpha))/lambda
  } else {
    ui <- rnorm(n,mean=0,sd=1)
  }
  ui <- ui*sigma + mu
  return(ui)
}

p.loggamma <- function(u,lambda) {
# derivarive wrt u of loggamma density
  csiLG(u,lambda)*dloggamma(x=u,lambda=lambda)
}

# Derivatives of cdf F_lambda(u) and density f_lambda(u)
# ======================================================

d.FL <- function(u, lambda) {
  zero <- 0.005
# first derivative with respect of lambda of the cdf
  if (abs(lambda) > zero) {
    temp <- try(integrate(intgpsiLG,lower=-50,upper=u,lambda=lambda), silent = TRUE)
    if (is.list(temp)) {
      res <- - temp$val
    } else {
      res <- NaN
    }
  } else { 
    if (lambda == 0) {
      res <- - integrate(intg3,lower=-6,upper=u)$val
    } else {
      ld <- -zero
      lo <-   0
      lu <-  zero 
      X  <- matrix(c(1,ld,ld^2,1,lo,lo^2,1,lu,lu^2),ncol=3,byrow=TRUE)
      y1 <- - integrate(intgpsiLG,lower=-50,upper=u,lambda=ld)$val
      y2 <- - integrate(intg3,lower=-6,upper=u)$val
      y3 <- - integrate(intgpsiLG,lower=-50,upper=u,lambda=lu)$val
      cf <- solve(X)%*%c(y1,y2,y3)
      res <- cf[1]+cf[2]*lambda+cf[3]*lambda^2
    }
  }
  return(res)
}

d2.FL <- function(u,lambda){ zero <- 0.005
# second derivative with respect of lambda of the cdf
if (abs(lambda) > zero) {res <- d2.FL.aux(u,lambda)}
else { 
if (lambda == 0) {i1 <- integrate(intg4,lower=-10,upper=u)$val
                  i2 <- integrate(intg6,lower=-10,upper=u)$val
                  res <- -i1 + i2}
else { y1 <- d2.FL.aux(u,-zero)
       y2 <- d2.FL.aux(u, zero)
       b   <- (y2-y1)/(2*zero)
       a   <- y1 - b*(-zero)
       res <- a+b*lambda  } }
res}

d.fL <- function(u,lambda){ zero=0.001
# first derivative with respect of lambda of the density function
if (abs(lambda) > zero)  {res <- d.fL.aux(u,lambda) }
if (abs(lambda) <= zero) {
       y1 <- d.fL.aux(u,-zero)
       y2 <- d.fL.aux(u, zero)
       b   <- (y2-y1)/(2*zero)
       a   <- y1 - b*(-zero)
       res <- a+b*lambda }
res}

# Auxiliary functions for derivatives

d2.FL.aux <- function(u,lambda){
i1 <- integrate(intgpsiLGd,lower=-50,upper=u,lambda=lambda)$val 
i2 <- integrate(intgpsiLG2,lower=-50,upper=u,lambda=lambda)$val 
-i1 + i2}

d.fL.aux <- function(u,lambda){
  lamm1  <-  lambda^(-1)
  lamm2  <-  lambda^(-2)
  lamm3  <-  lambda^(-3)
  llamm2 <-  log(lamm2)
  a1     <-  log(abs(lambda))+lamm2*llamm2+lamm1*u-lamm2*exp(lambda*u)
  a1d    <-  lamm1-2*lamm3*llamm2-2*lamm3-lamm2*u+2*lamm3*exp(lambda*u)-lamm2*u*exp(lambda*u)
  lres   <-  a1 - lgamma(lamm2) 
  res    <-  exp(lres)*(a1d+2*lamm3*digamma(lamm2))
res}

intg3 <- function(x){
(x^3/6)*dnorm(x,mean=0,sd=1)}

intgpsiLG <- function(x,lambda){
psiLG(x,lambda)*dloggamma(x=x, lambda=lambda)}


intgpsiLGd <- function(x,lambda){
psiLG.dot(x,lambda)*dloggamma(x=x, lambda=lambda)}

intgpsiLG2 <- function(x,lambda){
psiLG(x,lambda)^2*dloggamma(x=x, lambda=lambda)}

intg4 <- function(x){
((x^4+2)/12)*dnorm(x) }

intg6 <- function(x){
(x^3/6)^2*dnorm(x)}

# score functions csiLG and psiLG and their derivatives
# =================================================

csiLG <- function(u,lambda) {
# csiLG_lambda
zero <- 0.0001
if(abs(lambda) > zero) { res <- (1-exp(lambda*u))/lambda }
else { res<- -u }
res}

csiLG.dot  <- function(u,lambda) {
# derivative of csiLG_lambda wrt lambda
zero <- 0.0001
if(abs(lambda) > zero) {
 lamm2 <- lambda^(-2)
 res <- lamm2*(exp(lambda*u)*(1-lambda*u) - 1) 
} else { res <- -u^2/2 }
res}

csiLG.prime <- function(u,lambda) {
# derivative of csiLG_lambda wrt u
-exp(lambda*u)}

dersig.csiLG <- function(u,lambda,sigma,zero=0.0001){
# derivative wrt sigma of csiLG(u/sigma,lambda)
if(abs(lambda) > zero) { res <- u*exp(lambda*u)/sigma }
else { res <- u/sigma }
res}

psiLG <- function(u,lambda){
# psiLG_lambda
zero=0.001
if (abs(lambda) > zero) {res    <- psiLG.aux(u,lambda) } 
else  { 
  if (lambda == 0) {res <- (u^3)/6}
  else { 
  y1  <- psiLG.aux(u,-zero)
  y2  <- psiLG.aux(u, zero)
  b   <- (y2-y1)/(2*zero)
  a   <- y1 - b*(-zero)
  res <- a+b*lambda }}
res}

psiLG.prime <- function(u,lambda){
# derivative of psiLG_lambda wrt u
-csiLG.dot(u,lambda)}

psiLG.dot <- function(u,lambda){
# derivative of psiLG_lambda wrt lambda
zero <- 0.03
if (abs(lambda) > zero) {res <- psiLG.dot.aux(u,lambda) }
else  { 
  if (lambda == 0) { res <- (u^4+2)/12 }
    else  { 
    ld <- -zero
    lo <-   0
    lu <-  zero 
#    X  <- matrix(c(1,ld,ld^2,1,lo,lo^2,1,lu,lu^2),ncol=3,byrow=TRUE)
    y1 <- psiLG.dot.aux(u,ld)
    y2 <- (u^4+2)/12
    y3 <- psiLG.dot.aux(u,lu)
#    cf <- solve(X)%*%c(y1,y2,y3)
    cf1 <- y2
    cf2 <- ( lu*lu*(y1-y2)-ld*ld*(y3-y2)) / (ld*lu*lu-lu*ld*ld)
    cf3 <- ( lu*(y1-y2)-ld*(y3-y2) )/(ld*ld*lu-lu*lu*ld)
    res <- cf1+cf2*lambda+cf3*lambda^2}}
res}
                                      
psiLG.aux <- function(u,lambda){
lamu   <- u*lambda
elamu  <- exp(lamu)
lamm1  <- lambda^(-1)
lamm2  <- lambda^(-2)
llamm2 <- log(lamm2)
lamm3  <- lambda^(-3)
arg1d  <- lamm1-2*lamm3*llamm2-2*lamm3-lamm2*u+elamu*(2*lamm3-lamm2*u)
digl   <- digamma(lamm2)
-(arg1d+2*lamm3*digl) }

psiLG.dot.aux <- function(u,lambda){
l <- length(u)
g <- rep(1,l) 
lamm2  <-  lambda^(-2)
trigl  <-  trigamma(lamm2)
digl   <-  digamma(lamm2)
bl = -lambda^2 - 12*log(abs(lambda)) + 10 + 2*lambda*u + exp(lambda*u)*(-6+4*lambda*u-lambda^2*u^2)
ff = - bl + 6*digl + 4*trigl/lambda^2
g[ff> 0] = ff[ff >0]
rslt <- (ff > 0)* exp(log(g)-4*log(abs(lambda))) + (ff <=0) * ff*exp(- 4*log(abs(lambda)))
rslt}

######################################
# Score functions eta1, eta2, eta3, ... and derivatives for asymptotic covariance

eta1 <- function(u,lambda) csiLG(u,lambda) #fp.over.f.loggamma(u,lambda)
eta2 <- function(u,lambda) csiLG(u,lambda)*u + 1 #fp.over.f.loggamma(u,lambda)*u
eta3 <- function(u,lambda) psiLG(u,lambda) #l.score.loggamma(u,lambda)

eta1f <- function(u,lambda) eta1(u,lambda)*dloggamma(x=u, lambda=lambda)
eta2f <- function(u,lambda) eta2(u,lambda)*dloggamma(x=u, lambda=lambda)
eta3f <- function(u,lambda) eta3(u,lambda)*dloggamma(x=u, lambda=lambda)

eta11f <- function(u,lambda) (eta1(u,lambda))^2*dloggamma(x=u, lambda=lambda)
eta22f <- function(u,lambda) (eta2(u,lambda))^2*dloggamma(x=u, lambda=lambda)
eta33f <- function(u,lambda) (eta3(u,lambda))^2*dloggamma(x=u, lambda=lambda)
eta12f <- function(u,lambda) (eta1(u,lambda))*(eta2(u,lambda))*dloggamma(x=u, lambda=lambda)
eta13f <- function(u,lambda) (eta1(u,lambda))*(eta3(u,lambda))*dloggamma(x=u, lambda=lambda)
eta23f <- function(u,lambda) (eta2(u,lambda))*(eta3(u,lambda))*dloggamma(x=u, lambda=lambda)

##########################

eta1.mu     <- function(u,lambda) - csiLG.prime(u,lambda)
eta1.sigma  <- function(u,lambda) - csiLG.prime(u,lambda)*u
eta1.lambda <- function(u,lambda)   csiLG.dot(u,lambda)
eta2.mu     <- function(u,lambda) - csiLG.prime(u,lambda)*u -csiLG(u,lambda)

eta2.sigma  <- function(u,lambda) - csiLG.prime(u,lambda)*u^2 -csiLG(u,lambda)*u
eta2.lambda <- function(u,lambda)   csiLG.dot(u,lambda)*u   
eta3.mu     <- function(u,lambda) - psiLG.prime(u,lambda)
eta3.sigma  <- function(u,lambda) - psiLG.prime(u,lambda)*u
eta3.lambda <- function(u,lambda)   psiLG.dot(u,lambda)

#########################

eta1.muf     <- function(u,lambda) - (csiLG.prime(u,lambda))*dloggamma(x=u, lambda=lambda)
eta1.sigmaf  <- function(u,lambda) - (csiLG.prime(u,lambda)*u)*dloggamma(x=u, lambda=lambda)
eta1.lambdaf <- function(u,lambda) csiLG.dot(u,lambda)*dloggamma(x=u, lambda=lambda)
eta2.muf     <- function(u,lambda) (- csiLG.prime(u,lambda)*u -csiLG(u,lambda))*dloggamma(x=u, lambda=lambda)
eta2.sigmaf  <- function(u,lambda) (- csiLG.prime(u,lambda)*u^2 -csiLG(u,lambda)*u)*dloggamma(x=u, lambda=lambda)
eta2.lambdaf <- function(u,lambda) csiLG.dot(u,lambda)*u*dloggamma(x=u, lambda=lambda)
eta3.muf     <- function(u,lambda) - (psiLG.prime(u,lambda))*dloggamma(x=u, lambda=lambda)
eta3.sigmaf  <- function(u,lambda) - (psiLG.prime(u,lambda)*u)*dloggamma(x=u, lambda=lambda)
eta3.lambdaf <- function(u,lambda) psiLG.dot(u,lambda)*dloggamma(x=u, lambda=lambda)

###################################
