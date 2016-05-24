###
### margpimom.R
###

#fimomNeg works both for th of class vector and class matrix
#fpimomNeg, fppimomNeg work only for vector th
fimomNeg <- function(th, XtX, ytX, phi, tau) {
  p <- ncol(XtX)
  if (p==1) th <- matrix(th,ncol=1)
  if (is.vector(th)) {
    ans <- as.vector(.5*(mahalanobis(th,center=rep(0,p),cov=XtX,inverted=TRUE) - 2*ytX %*% matrix(th,ncol=1))/phi + tau*phi*sum(1/th^2) + sum(log(th^2)))
  } else {
    ans <- as.vector(.5*(mahalanobis(th,center=rep(0,p),cov=XtX,inverted=TRUE) - 2*ytX %*% t(th))/phi + tau*phi*rowSums(1/th^2) + rowSums(log(th^2)))
  }
  return(ans)
}

fpimomNeg <- function(th, XtX, ytX, phi, tau) (XtX %*% matrix(th,ncol=1) - t(ytX))/phi - 2*tau*phi/th^3 + 2/th
fppimomNeg <- function(th, XtX, ytX, phi, tau) XtX/phi + diag(6*tau*phi/th^4 - 2/th^2, ncol=length(th))

imomModeK <- function(thini, XtX, ytX, phi, tau) {
  #piMOM mode when phi is known using gradient algorithm
  th <- thini
  err <- 1; niter <- 0
  while ((err > 0.001) & (niter<50)) {
    err <- 0; niter <- niter+1
    for (i in 1:length(th)) {
      a <- c(2*tau*phi, 0, -2, (ytX[i]-sum(XtX[i,-i]*th[-i]))/phi, -XtX[i,i]/phi)
      thnew <- polyroot(a)
      thnew <- Re(thnew[abs(Im(thnew))< 1e-7])
      thnew <- thnew[sign(thnew)==sign(th[i])]
      err <- err+abs(th[i]-thnew)
      th[i] <- thnew
    }
  }
  return(th)  
}

imomIntegralApprox <- function(XtX, ytX, phi, tau, logscale=TRUE) {
#Laplace approx to product imom marginal (uses gradient search to find mode)
  m <- as.vector(solve(XtX + tau*diag(nrow(XtX))) %*% t(ytX))
  m <- imomModeK(m, XtX=XtX, ytX=ytX, phi=phi, tau=tau)
  V <- fppimomNeg(m, XtX=XtX, ytX=ytX, phi=phi, tau=tau)
  fopt <- fimomNeg(m,XtX=XtX,ytX=ytX,phi=phi,tau=tau)
  ans <- -fopt - .5*as.numeric(determinant(V,logarithm=TRUE)$modulus)
  if (!logscale) ans <- exp(ans)
  return(list(ans=ans,thopt=m,Vopt=V,objective=fopt))
}

#imomIntegralApproxOld <- function(XtX, ytX, phi, tau, logscale=TRUE) {
##Laplace approx to product imom marginal (uses numerical maximizer to find mode)
#  m <- as.vector(solve(XtX + tau*diag(nrow(XtX))) %*% t(ytX))
#  opt <- nlminb(m, objective=fimomNeg, XtX=XtX, ytX=ytX, phi=phi, tau=tau)
#  V <- fppimomNeg(opt$par, XtX=XtX, ytX=ytX, phi=phi, tau=tau)
#  ans <- -opt$objective - .5*as.numeric(determinant(V,logarithm=TRUE)$modulus)
#  if (!logscale) ans <- exp(ans)
#  return(list(ans=ans,thopt=opt$par,Vopt=V,objective=opt$objective))
#}


pimomMarginalK <- function(sel, y, x, phi, tau=1, method='Laplace', B=10^5, logscale=TRUE, XtX, ytX) {
#Marginal density of the data y~N(x*theta,phi*I) under a product imom prior (known variance)
# - sel: vector with indexes of variables included in the model
# - y: response variable
# - x: design matrix
# - phi: residual variance
# - tau: prior dispersion parameter
# - method: method to approximate the integral. 'Laplace' for Laplace approx. 'MC' to correct Laplace approx via Importance Sampling based on multivariate Cauchy
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
  phi <- as.double(phi); tau <- as.double(tau)
  if (method=='Laplace') method=0 else if (method=='MC') method=1 else if (method=='plugin') method=2 else stop("Invalid argument 'method'")
  method <- as.integer(method); B <- as.integer(B); logscale <- as.integer(logscale)
  ans <- .Call("pimomMarginalKI", sel, nsel, n, p, y, sumy2, XtX, ytX, phi, tau, method, B, logscale)
  return(ans)
}

pimomMarginalKR <- function(y, x, phi, tau=1, method='Laplace', B=10^5, logscale=TRUE, XtX, ytX) {
#Marginal density of the data y~N(x*theta,phi*I) under a product imom prior (known variance)
# - y: response variable
# - x: design matrix
# - phi: residual variance
# - tau: prior dispersion parameter
# - method: method to approximate the integral. 'Laplace' for Laplace approx. 'MC' to correct Laplace approx via Importance Sampling based on multivariate Cauchy
# - B: number of Monte Carlo samples to use (ignored unless method=='MC')
# - XtX, ytX: optionally, X'X and y'X can be specified to speed up computations
  if (is.matrix(y)) y <- as.vector(y)
  if (is.vector(x)) x <- matrix(x,ncol=1)
  p1 <- ncol(x); n <- nrow(x)
  if (p1==0) {
    ans <- sum(dnorm(y,0,sd=sqrt(phi),log=TRUE))
  } else {
    if (n != length(y)) stop("Dimensions of x and y don't match")
    if (missing(XtX)) { XtX <- t(x) %*% x }
    if (missing(ytX)) { ytX <- t(y) %*% x }
    ILaplace <- imomIntegralApprox(XtX=XtX,ytX=ytX,phi=phi,tau=tau,logscale=TRUE)
    k <- .5*p1*log(tau) - .5*sum(y^2)/phi - .5*n*log(2*pi) - .5*(n-p1)*log(phi) - .5*p1*log(pi)
    if (method=='Laplace') {
      ans <- k + ILaplace$ans
    } else if (method=='MC') {
      Vinv <- solve(ILaplace$Vopt)
      uplim <- ILaplace$thopt + 2*sign(ILaplace$thopt)*sqrt(diag(Vinv))
      sdprop <- diag(abs(uplim)/2,ncol=p1)
      Vprop <- sdprop %*% cov2cor(Vinv) %*% sdprop
      #thsim <- rmvnorm(B,rep(0,p1),Vprop)
      #adj <- - fimomNeg(thsim,XtX=XtX,ytX=ytX,phi=phi,tau=tau) - dmvnorm(thsim,rep(0,p1),Vprop,log=TRUE)
      thsim <- rmvt(B,sigma=Vprop,df=1) 
      adj <- - fimomNeg(thsim,XtX=XtX,ytX=ytX,phi=phi,tau=tau) - dmvt(thsim,delta=rep(0,p1),sigma=Vprop,df=1,log=TRUE)
      m <- max(adj)
      adj <- log(mean(exp(adj-m+500))) + m - 500
      ans <- k + adj
    } else {
      stop("Only method=='Laplace' and method=='MC' are currently implemented")
    }
  }
  if (!logscale) ans <- exp(ans)
  return(ans)
}



pimomMarginalU <- function(sel, y, x, alpha=0.001, lambda=0.001, tau=1, method='Laplace', B=10^5, logscale=TRUE, XtX, ytX) {
#Marginal density of the data y~N(x*theta,phi*I) under a product imom prior (unknown variance)
# - sel: vector with indexes of variables included in the model
# - y: response variable
# - x: design matrix
# - alpha, lambda: prior for phi is IGamma(alpha/2,lambda/2)
# - tau: prior dispersion parameter
# - method: method to approximate the integral for known phi. Integral wrt phi is performed via integrate. 'Laplace' for Laplace approx which may underestimate true value, 'MC' for exact evaluation which can be very computationally expensive. 'Hybrid' combines Laplace for fixed phi with numerical integration wrt phi. The Laplace error is estimated using a single exact evaluation for a value of phi close to the posterior mode
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
  tau <- as.double(tau)
  if (method=='Laplace') method=0 else if (method=='MC') method=1 else if (method=='plugin') method=2 else if (method=='Hybrid') method=3 else stop("Invalid argument 'method'")
  method <- as.integer(method); B <- as.integer(B); logscale <- as.integer(logscale)
  alpha <- as.double(alpha); lambda <- as.double(lambda)
  ans <- .Call("pimomMarginalUI",sel,nsel,n,p,y,sumy2,x,XtX,ytX,tau,method,B,logscale,alpha,lambda)
  return(ans);
}


fimomUNeg <- function(th, XtX, ytX, sumy2, tau, alpha, lambda, n) {
  #Last component in th is eta=log(phi) i.e. log-residual variance
  eta <- th[length(th)]; th <- th[-length(th)]; p <- length(th)
  ss <- lambda + sumy2 - 2*ytX %*% matrix(th,ncol=1) + mahalanobis(th,center=rep(0,p),cov=XtX,inverted=TRUE)
  ans <- tau*exp(eta)*sum(1/th^2) + sum(log(th^2)) + .5*eta*(n-p+alpha) + .5*exp(-eta)*ss
  return(ans)
}

fppimomUNeg <- function(th, XtX, ytX, sumy2, tau, alpha, lambda) {
  eta <- th[length(th)]; th <- th[-length(th)]; p <- length(th)
  ss <- lambda + sumy2 - 2*ytX %*% matrix(th,ncol=1) + mahalanobis(th,center=rep(0,p),cov=XtX,inverted=TRUE)
  ans <- matrix(NA,nrow=p+1,ncol=p+1)
  ans[1:p,1:p] <- exp(-eta) * XtX + diag(6*tau*exp(eta)/th^4 - 2/th^2,ncol=p)
  ans[1:p,p+1] <- ans[p+1,1:p] <- -exp(eta)*2*tau/th^3 - exp(-eta)*as.vector(matrix(th,nrow=1) %*% XtX - ytX)
  ans[p+1,p+1] <- tau*exp(eta)*sum(1/th^2) + .5*exp(-eta)*ss
  return(ans)
}

imomModeU <- function(thini, phiini, XtX, ytX, sumy2, tau, alpha, lambda, n) {
  #piMOM mode when phi is unknown using gradient algorithm
  th <- thini; phi <- phiini
  err <- 1; niter <- 0
  b <- -(n-ncol(XtX)+alpha)/2
  while ((err > 0.001) & (niter<50)) {
    err <- 0; niter <- niter+1
    for (i in 1:length(th)) {
      a <- c(2*tau*phi, 0, -2, (ytX[i]-sum(XtX[i,-i]*th[-i]))/phi, -XtX[i,i]/phi)
      thnew <- polyroot(a)
      thnew <- Re(thnew[abs(Im(thnew))< 1e-7])
      thnew <- thnew[sign(thnew)==sign(th[i])]
      err <- err+abs(th[i]-thnew)
      th[i] <- thnew
    }
    a <- tau*sum(1/th^2)
    b <- .5*(n-ncol(XtX)+alpha)
    c <- -.5*(lambda + sumy2 - 2*ytX %*% matrix(th,ncol=1) + matrix(th,nrow=1) %*% XtX %*% matrix(th,ncol=1))
    d <- sqrt(b^2 - 4*a*c)
    if (-b > d) { phinew <- (-b-d)/(2*a) } else { phinew <- (-b+d)/(2*a) }
    err <- err+abs(phi-phinew)
    phi <- phinew
  }
  return(c(th,log(phi)))  
}

imomUIntegralApprox <- function(thini, etaini, XtX, ytX, sumy2, tau, alpha, lambda, n, logscale=TRUE) {
#Laplace approx to product imom marginal (uses gradient search to find mode)
  m <- imomModeU(thini=thini,phiini=exp(etaini),XtX=XtX,ytX=ytX,sumy2=sumy2,tau=tau,alpha=alpha,lambda=lambda,n=n)
  fopt <- fimomUNeg(th=m,XtX=XtX,ytX=ytX,sumy2=sumy2,tau=tau,alpha=alpha,lambda=lambda,n=n)
  V <- fppimomUNeg(m,XtX=XtX,ytX=ytX,sumy2=sumy2,tau=tau,alpha=alpha,lambda=lambda)
  ans <- -fopt - .5*as.numeric(determinant(V,logarithm=TRUE)$modulus) + .5*nrow(XtX)*log(2*tau)
  if (!logscale) ans <- exp(ans)
  return(list(ans=ans,thopt=m,Vopt=V,objective=fopt))
}


imomUIntegralApproxOld <- function(thini, etaini, XtX, ytX, sumy2, tau, alpha, lambda, n, logscale=TRUE) {
#Laplace approx to product imom marginal (uses numerical optimizer to find mode)
  opt <- nlminb(c(thini,etaini), objective=fimomUNeg, XtX=XtX, ytX=ytX, sumy2=sumy2, tau=tau, alpha=alpha, lambda=lambda, n=n)
  V <- fppimomUNeg(opt$par,XtX=XtX,ytX=ytX,sumy2=sumy2,tau=tau,alpha=alpha,lambda=lambda)
  ans <- -opt$objective - .5*as.numeric(determinant(V,logarithm=TRUE)$modulus) + .5*nrow(XtX)*log(2*tau)
  if (!logscale) ans <- exp(ans)
  return(list(ans=ans,thopt=opt$par,Vopt=V,objective=opt$objective))
}


pimomMarginalUR <- function(y, x, alpha=0.001, lambda=0.001, tau=1, method='Laplace', B=10^5, logscale=TRUE, XtX, ytX) {
#Marginal density of the data y~N(x*theta,phi*I) under a product imom prior (unknown variance)
# - y: response variable
# - x: design matrix
# - alpha, lambda: prior for phi is IGamma(alpha/2,lambda/2)
# - tau: prior dispersion parameter
# - method: method to approximate the integral for known phi. Integral wrt phi is performed via integrate. 'Laplace' for Laplace approx which may underestimate true value, 'MC' for exact evaluation which can be very computationally expensive. 'Hybrid' uses numerical integration to integrate over phi and Laplace to integrate over theta (it corrects the Laplace error with an exact evaluation for a single value of phi close to the posterior mode)
# - B: number of Monte Carlo samples to use (ignored if method=='Laplace')
  #require(actuar)
  if (is.matrix(y)) y <- as.vector(y)
  if (is.vector(x)) x <- matrix(x,ncol=1)
  if (ncol(x)==0) {
    n <- length(y)
    term1 <- .5*(n + alpha)
    num <- .5*alpha*log(lambda) + lgamma(term1)
    den <- .5*n*log(pi) + lgamma(alpha/2)
    ans <- num -den - term1*log(lambda + sum(y^2))
  } else {
    if (missing(XtX)) { XtX <- t(x) %*% x }
    if (missing(ytX)) { ytX <- t(y) %*% x }
    if (method=='Hybrid') {
      f2int <- function(z,method, adj=1) {
        ans <- double(length(z))
        for (i in 1:length(ans)) ans[i] <- pimomMarginalKR(y=y,x=x,phi=z[i],tau=tau,method=method,B=B,logscale=FALSE,XtX=XtX,ytX=ytX)
        ans <- ans * dinvgamma(z, alpha/2, scale=lambda/2) * adj
        return(ans)
      }
      #Compute adjustment factor
      e <- y - x %*% solve(XtX + diag(tau,nrow=nrow(XtX))) %*% t(ytX)
      phiest <- (sum(e^2)+lambda)/(length(y)+alpha)
      intmc <- f2int(phiest,method='MC')
      intlapl <- f2int(phiest,method='Laplace')
      adj <- intmc/intlapl
      ans <- log(integrate(f2int, 0, Inf, method='Laplace', adj=adj)$value)
    } else {
      thini <- as.vector(solve(XtX + tau*diag(nrow(XtX))) %*% t(ytX))
      e <- y - x %*% matrix(thini,ncol=1)
      etaini <- log((sum(e^2)+lambda)/(length(y)+alpha))
      ans <- imomUIntegralApprox(thini=thini,etaini=etaini,XtX=XtX,ytX=ytX,sumy2=sum(y^2),tau=tau,alpha=alpha,lambda=lambda,n=nrow(x),logscale=TRUE)
      ans <- ans$ans + .5*alpha*log(.5*lambda) - .5*nrow(x)*log(2*pi) - lgamma(.5*alpha)
    }
  }
  if (!logscale) ans <- exp(ans)
  return(ans)
}

