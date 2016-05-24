###
### imomunknown.R
###

imomunknown <- function(theta1hat,V1,n,nuisance.theta,g=1,nu=1,theta0,ssr,method='adapt',nquant=100,B=10^5) {

fncp <- function(sigma2) {
  l <- theta1hat-theta0; l <- matrix(l,nrow=1) %*% solve(V1) %*% matrix(l,ncol=1) / sigma2
  l[l==Inf] <- exp(80)
  return(l)
}
f <- function(z) { ans <- (n*gi/z)^((nu+p1)/2) * exp(-n*gi/z); ans[z==0] <- 0; return(ans) }

if (missing(theta0)) theta0 <- rep(0,length(theta1hat))
p1 <- length(theta1hat)
m <- double(length(g))
if (method=='MC') {
  sigma2 <- 1/rgamma(B,(n-nuisance.theta)/2,ssr/2)
  l <- fncp(sigma2)
  z <- rchisq(B,df=p1,ncp=l)
  for (i in 1:length(m)) { gi <- g[i]; m[i] <- mean(f(z)) }
} else if (method=='adapt') {
  f2 <- function(z,z2) { return(f(z)*dchisq(z,df=p1,ncp=fncp(z2))) }
  qseq <- 1/qgamma((2*(1:nquant)-1)/(2*nquant),(n-nuisance.theta)/2,ssr/2)
  for (i in 1:length(m)) {
    m[i] <- 0; for (j in 1:nquant) { gi <- g[i]; m[i] <- m[i]+integrate(f2,0,Inf,z2=qseq[j])$value }
    m[i] <- m[i]/nquant
  }
} else {
  stop('method must be adapt or MC')
}
t1 <- theta1hat-theta0; t1 <- matrix(t1,nrow=1) %*% solve(V1) %*% matrix(t1,ncol=1) / ssr #noncentr param
bf <- exp((p1/2)*log(2/(n*g)) + lgamma(p1/2) - lgamma(nu/2) + ((n-nuisance.theta)/2)*log(1+t1)) * m
return(bf)
}

