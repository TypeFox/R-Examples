###
### imombf.lm.R
###

imombf.lm <- function(lm1,coef,g,prior.mode,nu=1,theta0,method='adapt',nquant=100,B=10^5) {
if ((!missing(g)) & (!missing(prior.mode))) warning('Both g and prior.mode were specified. g will be ignored')
if ((missing(g)) & (missing(prior.mode))) stop('Either g or prior.mode must be specified')
if (missing(theta0)) theta0 <- rep(0,length(coef)) else if (length(theta0)!=length(coef)) stop('theta0 must have the same length as coef')
  
  thetahat <- coef(lm1)
  V <- summary(lm1)$cov.unscaled
  n <- length(lm1$residuals); p <- length(thetahat); p1 <- length(coef)
  if ((min(coef)<1) | (max(coef)>p)) stop('Non-valid value for coef. Use only values between 1 and the number of coefficients in lm1')
  ssr <- sum(residuals(lm1)^2); sr <- sqrt(ssr/(n-p))
  if (!missing(prior.mode)) g <- mode2g(prior.mode,prior='iMom')
  bf.imom <- imomunknown(thetahat[coef],V[coef,coef],n,nuisance.theta=p-p1,g=g,nu=nu,theta0=theta0,ssr=ssr,method=method,nquant=nquant,B=B)
  return(bf.imom)
}

