###
### margpmom.R
###

eprod <- function(m, S, power=1, dof= -1) {
  #Mean of prod (x_i)^(2*power) when x_i ~ T_dof(mu,sigma). Set dof=-1 for N(mu,sigma). Written by John Cook
  ans <- .Call("eprod_I",as.double(m),as.double(S), as.integer(length(m)), as.integer(power), as.double(dof))
  ans
}

fmomNeg <- function(th, m, S, phi, tau, r, logscale=TRUE) .5*mahalanobis(th, center=m, cov=S, inverted=TRUE)/phi - r*sum(log(th^2))
fpmomNeg <- function(th, m, S, phi, tau, r) S %*% matrix(th-m, ncol=1)/phi - 2*r/th
fppmomNeg <- function(th, m, S, phi, tau, r) S/phi + 2*r*diag(1/th^2,nrow=length(th))

pmomIntegralApproxR <- function(m, S, phi, tau, r, logscale=TRUE) {
  #Laplace approx to integral N(th; m, phi*solve(S)) prod (th/(phi*tau))^2r wrt th
  opt <- nlminb(m, objective=fmomNeg, gradient=fpmomNeg, m=m, S=S, phi=phi, tau=tau, r=r)$par
  fopt <- -fmomNeg(opt,m=m,S=S,phi=phi,tau=tau,r=r)
  hess <- fppmomNeg(opt,m=m,S=S,phi=phi,tau=tau,r=r)
  ans <- fopt + .5*log(det(S)) - .5*log(det(hess)) - .5*length(m)*log(phi)
  if (!logscale) ans <- exp(ans)
  return(ans)
}

pmomMarginalK <- function(sel, y, x, phi, tau, r=1, method='auto', B=10^5, logscale=TRUE, XtX, ytX) {
#Marginal density of the data y~N(x*theta,phi*I) under a product mom prior (known variance)
# - sel: vector with indexes of variables included in the model
# - y: response variable
# - x: design matrix
# - phi: residual variance
# - tau: prior dispersion parameter
# - r: prior power parameter is 2*r
# - method: method to approximate the integral. 'Laplace' for Laplace approx. 'MC' for Monte Carlo. 'Plug-in' for plug-in estimate. 'auto' for exact calculation if p<=10, else Laplace approx
# - B: number of Monte Carlo samples to use (ignored unless method=='MC')
# - XtX, ytX: optionally, X'X and y'X can be specified to speed up computations
  if (is.matrix(y)) y <- as.vector(y)
  if (is.vector(x)) x <- matrix(x,ncol=1)
  if (missing(XtX)) { XtX <- t(x) %*% x } else { XtX <- as.matrix(XtX) }
  if (missing(ytX)) { ytX <- as.vector(matrix(y,nrow=1) %*% x) } else { ytX <- as.vector(ytX) }
  if (is.logical(sel)) sel <- which(sel)
  if ((length(sel)>0) && ((min(sel)<1) | max(sel)>ncol(x))) stop('Invalid specification of parameter sel')
  sel <- as.integer(sel-1); nsel <- as.integer(length(sel)); 
  p <- as.integer(ncol(x)); n <- as.integer(nrow(x))
  y <- as.double(y); sumy2 <- sum(y^2)
  phi <- as.double(phi); tau <- as.double(tau); r <- as.integer(r)
  if (method=='auto') method=-1 else if (method=='Laplace') method=0 else if (method=='MC') method=1 else if (method=='plugin') method=2 else stop("Invalid 'method'")
  method <- as.integer(method)
  B <- as.integer(B); logscale <- as.integer(logscale)
  ans <- .Call("pmomMarginalKI", sel, nsel, n, p, y, sumy2, XtX, ytX, phi, tau, r, method, B, logscale)
  return(ans)
}


pmomMarginalKR <- function(y, x, phi, tau, r=1, method='Laplace', B=10^5, logscale=TRUE) {
  #Marginal likelihood for product moment prior (variance phi known)
  # - Likelihood: y ~ N(x %*% th, phi * I)
  # - Prior proportional to N(th; 0, tau*phi*I) * prod(th^2/(phi*tau))^r
  #   i.e. phi is the residual variance; tau the prior dispersion parameter
  #require(mvtnorm)
  n <- length(y); p <- ncol(x)
  if (p==0) {
    ans <- sum(dnorm(y,0,sd=sqrt(phi),log=TRUE))
  } else {
    S <- t(x) %*% x + diag(p)/tau
    m <- solve(S) %*% t(x) %*% matrix(y,ncol=1)
    ans <- -.5*(sum(y^2) - t(m) %*% S %*% m)/phi - p*(sum(log(seq(1,2*r-1,by=2)))) - .5*n*log(2*pi*phi) - .5*p*log(tau) - log(sqrt(det(S))) - r*p*log(tau*phi)
    if (method=='Laplace') {
      I <- pmomIntegralApproxR(m=m, S=S, phi=phi, tau=tau, r=r, logscale=TRUE)
    } else if (method=='1storder') {
      I <- r*sum(log(m^2))
    } else if (method=='MC') {
      thsim <- rmvnorm(B,m,phi*solve(S))
      eprod <- exp(rowSums(log(thsim^(2*r))))
      I <- log(mean(eprod))
    }
    ans <- ans + I
  }
  if (!logscale) ans <- exp(ans)
  return(ans)
}


pmomMarginalU <- function(sel, y, x, alpha=0.001, lambda=0.001, tau=1, r=1, method='auto', B=10^5, logscale=TRUE, XtX, ytX) {
#Marginal density of the data y~N(x*theta,phi*I) under a product imom prior (unknown variance)
# - sel: vector with indexes of variables included in the model
# - y: response variable
# - x: design matrix
# - alpha, lambda: prior for phi is IGamma(alpha/2,lambda/2)
# - tau: prior dispersion parameter
# - r: prior power parameter is 2*r
# - method: method to approximate the integral. 'Laplace' for Laplace approx. 'MC' for Monte Carlo. 'Plug-in' for plug-in estimate.
# - B: number of Monte Carlo samples to use (ignored unless method=='MC')
  if (is.matrix(y)) y <- as.vector(y)
  if (is.vector(x)) { x <- matrix(x,ncol=1) } else { x <- as.matrix(x) }
  if (missing(XtX)) { XtX <- t(x) %*% x } else { XtX <- as.matrix(XtX) }
  if (missing(ytX)) { ytX <- as.vector(matrix(y,nrow=1) %*% x) } else { ytX <- as.vector(ytX) }
  if (is.logical(sel)) sel <- which(sel)
  if ((length(sel)>0) && ((min(sel)<1) | max(sel)>ncol(x))) stop('Invalid specification of parameter sel')
  sel <- as.integer(sel-1); nsel <- as.integer(length(sel)); 
  p <- as.integer(ncol(x)); n <- as.integer(nrow(x))
  y <- as.double(y); sumy2 <- sum(y^2)
  tau <- as.double(tau); r <- as.integer(r)
  if (method=='auto') method=-1 else if (method=='Laplace') method=0 else if (method=='MC') method=1 else if (method=='plugin') method=2 else stop("Invalid 'method'")
  method <- as.integer(method)
  B <- as.integer(B); logscale <- as.integer(logscale)
  alpha <- as.double(alpha); lambda <- as.double(lambda)
  ans <- .Call("pmomMarginalUI",sel,nsel,n,p,y,sumy2,x,XtX,ytX,tau,r,method,B,logscale,alpha,lambda)
  return(ans);
}

pmomMarginalUR <- function(y, x, r, alpha=0.001, lambda=0.001, tau, method='Laplace', B=10^5, logscale=TRUE) {
  #Marginal likelihood for product moment prior (variance phi unknown)
  # - Likelihood: y ~ N(x %*% th, phi * I)
  # - Prior for th proportional to N(th; 0, tau*phi*I) * prod(th^2/(phi*tau))^r
  # - Prior for phi: IGamma(alpha/2,lambda/2)
  #   i.e. phi is the residual variance; tau the prior dispersion parameter
  #require(mvtnorm)
  n <- length(y); p <- ncol(x)
  if (ncol(x)==0) {
    term1 <- .5*(n + alpha)
    num <- .5*alpha*log(lambda) + lgamma(term1)
    den <- .5*n*log(pi) + lgamma(alpha/2)
    ans <- num -den - term1*log(lambda + sum(y^2))
  } else {
    S <- t(x) %*% x + diag(p)/tau
    m <- solve(S) %*% t(x) %*% matrix(y,ncol=1)
    nu <- 2*r*p + n + alpha
    ss <- as.numeric(lambda + sum(y^2) - t(m) %*% S %*% m)
    V <- S*nu/ss
    #
    if (method=='Laplace') {
      I <- pmomIntegralApproxR(m=m, S=S, phi=nu/(nu-2), tau=tau, r=r, logscale=TRUE)
    } else if (method=='1storder') {
      I <- r*sum(log(m^2))
    } else if (method=='MC') {
      cholV <- t(chol(solve(V)))
      z <- rmvnorm(B,rep(0,p),diag(p))
      thsim <- as.vector(m) + (cholV %*% t(z)) * sqrt(nu/rchisq(B,df=nu))
      eprod <- exp(colSums(log(thsim^(2*r))))
      seprod <- sd(eprod)/sqrt(length(eprod))
      I <- log(mean(eprod))
    } else {
      stop("Only 'Laplace', '1storder' and 'MC' methods are implemented")
    }
    #
    num <- lgamma(nu/2) + .5*alpha*log(lambda/2) + .5*nu*log(2) - .5*nu*log(ss)
    den <- p*(sum(log(seq(1,2*r-1,by=2)))) + .5*n*log(2*pi) + .5*log(det(S)) + (.5*p+r*p)*log(tau) + lgamma(alpha/2)
    ans <- I + num - den
  }
  if (!logscale) ans <- exp(ans)
  return(ans)
}
