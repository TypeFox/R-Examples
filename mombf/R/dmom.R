###
### dmom.R
###

#Wrapper to call dpmom (product mom) and dqmom (quadratic mom)
dmom <- function(x, tau, a.tau, b.tau, phi=1, r=1, V1, baseDensity='normal', nu=3, logscale=FALSE, penalty='product') {
  if (penalty=='product') {
    ans <- dpmom(x, tau=tau, a.tau=a.tau, b.tau=b.tau, phi=phi, r=r, baseDensity=baseDensity, logscale=logscale)
  } else if (penalty=='quadratic') {
    if (r>1) stop("r>1 not implemented for penalty=='quadratic'. Try penalty=='product' instead")
    ans <- dqmom(x, V1=V1, g=tau, n=1, baseDensity=baseDensity, nu=nu)
  } else {
    stop("Only 'penalty==product' and 'penalty==quadratic' are implemented")
  }
  return(ans)
}

##Product MOM density

setMethod("dpmom", signature(x='vector'), function(x, tau, a.tau, b.tau, phi=1, r=1, baseDensity='normal', logscale=FALSE) {
  if (baseDensity!='normal') stop("Only baseDensity=='normal' is implemented for the product MOM")
  if (missing(tau) & (missing(a.tau) | missing(b.tau))) stop("Either tau or (a.tau,b.tau) must be specified")
  if (!missing(tau)) {
    ans <- dnorm(x,0,sd=sqrt(tau*phi),log=TRUE) + r*log(x^2/(tau*phi)) - sum(log(seq(1,2*r-1,by=2)))
  } else {
    ct <- lgamma(r+.5*(a.tau+1)) + r*log(2) - lgamma(.5*a.tau) - .5*log(pi) - (r+.5)*log(b.tau) - .5*log(phi)
    ans <- ct + r*log(x^2/phi) - sum(log(seq(1,2*r-1,by=2))) - (r+.5*(a.tau+1))*log(1+x^2/(b.tau*phi))
  }
  if (!logscale) ans <- exp(ans)
  ans
}
)
setMethod("dpmom", signature(x='matrix'), function(x, tau, a.tau, b.tau, phi=1, r=1, baseDensity='normal', logscale=FALSE) {
  if (baseDensity!='normal') stop("Only baseDensity=='normal' is implemented for the product MOM")
  if (missing(tau) & (missing(a.tau) | missing(b.tau))) stop("Either tau or (a.tau,b.tau) must be specified")
  p <- ncol(x)
  normct <- p*sum(log(seq(1,2*r-1,by=2)))
  distval <- rowSums(x^2)
  if (!missing(tau)) {
    ans <- -(p * log(2 * pi) + p*(log(phi)+log(tau))  + distval/(phi*tau))/2 + r*rowSums(log(x^2/(tau*phi))) - normct
  } else {
    anew <- r*p + .5*p + .5*a.tau
    num <- lgamma(anew) - anew*log(1+rowSums(x^2)/(phi*b.tau))  + r*rowSums(log(2*x^2/(phi*b.tau))) - normct
    den <- .5*p*(log(pi)+log(phi)+log(b.tau)) + lgamma(.5*a.tau)
    ans <- num - den
  }
  if (!logscale) ans <- exp(ans)
  ans
}
)

##Quadratic MOM density
dqmom <- function(x,V1=1,g,n=1,theta0,baseDensity='normal',nu=3) {
if (missing(g)) stop("Prior dispersion must be specified for quadratic MOM")
if (!(baseDensity %in% c('normal','t'))) stop("The only implemented baseDensity values are 'normal' and 't'")
if (baseDensity=='t' & nu<3) stop('nu must be >=3, otherwise the prior is improper')
if (missing(V1)) {
  if (is.vector(x)) V1 <- 1 else V1 <- diag(ncol(x))
}
if (missing(theta0)) {
  if (is.vector(V1)) theta0 <- 0 else theta0 <- rep(0,ncol(V1))
}
    
if (is.vector(V1)) {
  qtheta <- (x-theta0)^2/(n*g*V1)
  if (baseDensity=='normal') {
      ans <- qtheta*dnorm(x,theta0,sd=sqrt(n*g*V1))
  } else if (baseDensity=='t') {
    normct <- exp(lgamma(.5*(nu+1))-lgamma(.5*nu)-.5*log(nu*pi*n*g*V1))
    ans <- qtheta*normct*(1+qtheta/nu)^(-.5*(nu+1))*(nu-2)/nu
  }
} else {
  #require(mvtnorm)
  qtheta <- mahalanobis(x,center=theta0,cov=n*g*V1)
  if (baseDensity=='normal') {
    ans <- qtheta*dmvnorm(x,mean=theta0,sigma=n*g*V1)/ncol(V1)
  } else if (baseDensity=='t') {
    ans <- qtheta*dmvt(x,delta=theta0,sigma=n*g*V1)*(nu-2)/(nu*length(theta0))
  }
}
return(ans)
}

