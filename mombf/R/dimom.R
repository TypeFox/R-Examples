###
### dimom.R
###

#Wrapper to call dpimom (product iMOM) and dqimom (quadratic imom)
dimom <- function(x, tau=1, phi=1, V1, logscale=FALSE, penalty='product') {
  if (penalty=='product') {
    ans <- dpimom(x, tau=tau, phi=phi, logscale=logscale)
  } else if (penalty=='quadratic') {
    ans <- dqimom(x, V1=V1, g=tau, n=1, nu=1, logscale=logscale)
  } else {
    stop("Only 'penalty==product' and 'penalty==quadratic' are implemented")
  }
  return(ans)
}

#Product iMOM

setMethod("dpimom", signature(x='vector'), function(x, tau=1, phi=1, logscale=FALSE) {
  ans <- .5*(log(tau)+log(phi)) - lgamma(.5) - log(x^2) - tau*phi/x^2
  ans[is.nan(ans)] <- -Inf
  if (!logscale) ans <- exp(ans)
  ans
}
)

setMethod("dpimom", signature(x='matrix'), function(x, tau=1, phi=1, logscale=FALSE) {
  x2 <- x^2
  ans <- ncol(x)*(.5*(log(tau)+log(phi)) - lgamma(.5)) - rowSums(log(x2)) - tau*phi*rowSums(1/x2)
  ans[is.nan(ans)] <- -Inf
  if (!logscale) ans <- exp(ans)
  ans
}
)

setMethod("dpimom", signature(x='data.frame'), function(x, tau=1, phi=1, logscale=FALSE) {
  dpimom(as.matrix(x),tau=tau,phi=phi,logscale=logscale)
}
)


#Quadratic iMOM
dqimom <- function(x,V1=1,g=1,n=1,nu=1,theta0,logscale=FALSE) {
if (is.vector(x)) {
  if (missing(theta0)) theta0 <- 0
  if (missing(V1)) V1 <- 1
  qtheta <- (x-theta0)^2/(n*g*V1)
  k <- -0.5*log(n*g*V1) - lgamma(nu/2)
  p1 <- 1
} else {
  if (missing(theta0)) theta0 <- rep(0,ncol(x))
  if (missing(V1)) V1 <- diag(ncol(x))
  qtheta <- t(matrix(x,nrow=nrow(x))) - theta0
  qtheta <- qtheta %*% solve(n*g*V1) %*% t(qtheta)
  k <- -0.5*log(n*g*det(V1)) - lgamma(nu/2) + lgamma(ncol(x)/2) - .5*ncol(x)*log(pi)
  p1 <- ncol(x)
}
ans <- (k - .5*(nu+p1)*log(qtheta) -1/qtheta)
if (!logscale) { ans <- exp(ans); ans[is.na(ans)] <- 0 }
return(ans)
}
