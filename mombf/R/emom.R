###
### emom.R
###


setMethod("demom",signature(x='vector'),function(x, tau, a.tau, b.tau, phi=1, logscale=FALSE) {
  V1 <- 1
  if (!missing(tau)) {
    pen <- -tau*phi/x^2
    normct <- sqrt(2)
    ans <- pen + dnorm(x,mean=0,sd=sqrt(tau*phi*V1),log=TRUE) + normct
  } else {
    p <- 1; x2phi <- x^2/phi
    anew <- .5*(a.tau+p); bnew <- .5*(b.tau+x2phi)
    num <- sqrt(2)*p + .5*a.tau*log(.5*b.tau)
    den <- lgamma(.5*a.tau) + .5*p*(log(2*pi)+log(phi)) + anew*log(bnew)
    bt <- bnew/x2phi
    ans <- num - den + log(2) + .5*anew*log(bt) + log(besselK(sqrt(4*bt),nu=anew,expon.scaled=TRUE)) - sqrt(4*bt)
  }
  if (!logscale) ans <- exp(ans)
  return(ans)
}
)

setMethod("demom",signature(x='data.frame'),function(x, tau, a.tau, b.tau, phi=1, logscale=FALSE) {
  demom(as.matrix(x),tau=tau,a.tau=a.tau,b.tau=b.tau,phi=phi,logscale=logscale)
}
)

setMethod("demom",signature(x='matrix'),function(x, tau, a.tau, b.tau, phi=1, logscale=FALSE) {
  p <- ncol(x)
  V1 <- diag(p)
  if (!missing(tau)) {
    pen <- -tau*phi*rowSums(1/x^2)
    normct <- p*sqrt(2)
    ans <- pen + dmvnorm(x,mean=rep(0,p),sigma=tau*phi*V1,log=TRUE) + normct
  } else {
    anew <- .5*(a.tau+p); bnew <- .5*(b.tau+rowSums(x^2)/phi)
    num <- sqrt(2)*p + .5*a.tau*log(.5*b.tau)
    den <- lgamma(.5*a.tau) + .5*p*(log(2*pi)+log(phi)) + anew*log(bnew)
    bt <- bnew*phi*rowSums(1/x^2)
    ans <- num - den + log(2) + .5*anew*log(bt) + log(besselK(sqrt(4*bt),nu=anew,expon.scaled=TRUE)) - sqrt(4*bt)
  }
  if (!logscale) ans <- exp(ans)
  return(ans)
}
)


pemom <- function(q, tau, a.tau, b.tau) integrate(demom,-Inf,q,tau=tau,a.tau=a.tau,b.tau=b.tau)$value

