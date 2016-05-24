###
### margpemom.R
###

femomNeg <- function(th, m, S, phi, tau, logscale=TRUE) .5*mahalanobis(th, center=m, cov=S, inverted=TRUE)/phi + tau*phi*sum(1/th^2)
fpemomNeg <- function(th, m, S, phi, tau) S %*% matrix(th-m, ncol=1)/phi - 2*tau*phi*sum(1/th^3)
fppemomNeg <- function(th, m, S, phi, tau) S/phi + 6*tau*phi*diag(1/th^4,nrow=length(th))

pemomIntegralApproxR <- function(m, S, phi, tau, logscale=TRUE) {
  #Laplace approx to integral N(th; m, phi*solve(S)) prod(exp(-tau*phi/th^2)) wrt th
  opt <- nlminb(m, objective=femomNeg, gradient=fpemomNeg, m=m, S=S, phi=phi, tau=tau)$par
  fopt <- -femomNeg(opt,m=m,S=S,phi=phi,tau=tau)
  hess <- fppemomNeg(opt,m=m,S=S,phi=phi,tau=tau)
  ans <- fopt + .5*log(det(S)) - .5*log(det(hess)) - .5*length(m)*log(phi)
  if (!logscale) ans <- exp(ans)
  return(ans)
}



pemomMarginalKR <- function(y, x, phi, tau, method='Laplace', B=10^5, logscale=TRUE) {
  #Marginal likelihood for product eMOM prior (variance phi known)
  # - Likelihood: y ~ N(x %*% th, phi * I)
  # - Prior proportional to N(th; 0, tau*phi*I) * prod(exp(-tau*phi/th^2)^r
  #   i.e. phi is the residual variance; tau the prior dispersion parameter
  #require(mvtnorm)
  n <- length(y); p <- ncol(x)
  if (p==0) {
    ans <- sum(dnorm(y,0,sd=sqrt(phi),log=TRUE))
  } else {
    S <- t(x) %*% x + diag(p)/tau
    m <- solve(S) %*% t(x) %*% matrix(y,ncol=1)
    ans <- -.5*(sum(y^2) - t(m) %*% S %*% m)/phi + p*sqrt(2)  - .5*n*log(2*pi*phi) - .5*p*log(tau) - log(sqrt(det(S))) 
    if (method=='Laplace') {
      I <- pemomIntegralApproxR(m=m, S=S, phi=phi, tau=tau, logscale=TRUE)
    } else if (method=='1storder') {
      I <- - tau*phi*sum(1/m^2)
    } else if (method=='MC') {
      thsim <- rmvnorm(B,m,phi*solve(S))
      eprod <- exp(-tau*phi*rowSums(1/thsim^2))
      I <- log(mean(eprod))
    }
    ans <- ans + I
  }
  if (!logscale) ans <- exp(ans)
  return(ans)
}




pemomMarginalUR <- function(y, x, r, alpha=0.001, lambda=0.001, tau, method='Laplace', B=10^5, logscale=TRUE) {
  #Marginal likelihood for product eMOM prior (variance phi unknown)
  # - Likelihood: y ~ N(x %*% th, phi * I)
  # - Prior for th proportional to N(th; 0, tau*phi*I) * prod(exp(-tau*phi/th^2)^r
  # - Prior for phi: IGamma(alpha/2,lambda/2)
  #   i.e. phi is the residual variance; tau the prior dispersion parameter
  #require(mvtnorm) 
  if (is.vector(x)) x <- matrix(x,ncol=1)
  n <- length(y); p <- ncol(x)
  if (ncol(x)==0) {
    term1 <- .5*(n + alpha)
    num <- .5*alpha*log(lambda) + lgamma(term1)
    den <- .5*n*log(pi) + lgamma(alpha/2)
    ans <- num -den - term1*log(lambda + sum(y^2))
  } else {
    S <- t(x) %*% x + diag(p)/tau
    m <- solve(S) %*% t(x) %*% matrix(y,ncol=1)
    apost <- n + alpha
    lpost <- as.numeric(lambda + sum(y^2) - t(m) %*% S %*% m)
    #
    if (method=='Laplace') {
      pen <- lpost*tau*sum(1/m^2)
      I <- log(2) - lgamma(.5*apost) + .25*apost*log(.5*pen) + lbesselK(sqrt(2*pen),nu=.5*apost)
    } else if (method=='1storder') {
      phi <- lpost/apost
      I <- - tau*phi*sum(1/m^2)
    } else if (method=='MC') {
      phisim <- 1/rgamma(B, apost, lpost)
      cholV <- t(chol(solve(S)))
      z <- rmvnorm(B,rep(0,p),diag(p))
      thsim <- as.vector(m) + (cholV %*% t(z)) * sqrt(phisim)
      eprod <- -tau*phisim*colSums(1/thsim^2)
      offset <- max(eprod)
      I <- log(mean(exp(eprod-offset))) + offset
    } else {
      stop("Only 'Laplace', '1storder' and 'MC' methods are implemented")
    }
    #
    num <- p*sqrt(2) + .5*alpha*log(lambda/2) + lgamma(apost/2)
    den <- .5*n*log(2*pi) + .5*log(det(S)) + .5*p*log(tau) + lgamma(alpha/2) + .5*apost*log(lpost/2)
    ans <- I + num - den
  }
  if (!logscale) ans <- exp(ans)
  return(ans)
}
