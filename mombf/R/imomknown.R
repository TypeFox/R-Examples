###
### imomknown.R.R
###

imomknown <- function(theta1hat,V1,n,nuisance.theta,g=1,nu=1,theta0,sigma,method='adapt',B=10^5) {
if (missing(sigma)) stop('sigma must be specified')
if (missing(theta0)) theta0 <- rep(0,length(theta1hat))
f <- function(z) { ans <- (n*gi/z)^((nu+p1)/2) * exp(-n*gi/z); ans[z==0] <- 0; return(ans) }

p1 <- length(theta1hat)
l <- theta1hat-theta0; l <- matrix(l,nrow=1) %*% solve(V1) %*% matrix(l,ncol=1) / sigma^2 #noncentr param
m <- double(length(g))
if (method=='MC') {
  z <- rchisq(B,df=p1,ncp=l)
  for (i in 1:length(m)) { gi <- g[i]; m[i] <- mean(f(z)) }
} else if (method=='adapt') {
  f2 <- function(z) { return(f(z)*dchisq(z,df=p1,ncp=l)) }
  for (i in 1:length(m)) { gi <- g[i]; m[i] <- integrate(f2,0,Inf)$value }
} else {
  stop('method must be adapt or MC')
}
bf <- exp((p1/2)*log(2/(n*g)) + lgamma(p1/2)-lgamma(nu/2) + .5*l) * m
return(bf)
}

