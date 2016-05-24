require(sde)

# ex3.08.R
set.seed(123)
d <- expression(-1*x)
s <- expression(2) 
sde.sim(drift=d, sigma=s,N=50,delta=0.01) -> X
S <- function(t, x, theta) sqrt(theta[2])
B <- function(t, x, theta) -theta[1]*x
 
true.loglik <- function(theta) {
 DELTA <- deltat(X)
 lik <- 0
 for(i in 2:length(X))
  lik <- lik + dnorm(X[i], mean=X[i-1]*exp(-theta[1]*DELTA), 
  sd = sqrt((1-exp(-2*theta[1]*DELTA))*theta[2]/(2*theta[1])),TRUE)
  lik  
}
 
xx <- seq(-10,10,length=20)
sapply(xx, function(x) true.loglik(c(x,4))) -> py
sapply(xx, function(x) EULERloglik(X,c(x,4),B,S)) -> pz
sapply(xx, function(x) SIMloglik(X,c(x,4),B,S,M=10000,N=5)) -> pw

par(mar=c(5,5,1,1))
plot(xx,py, type="l", xlab=expression(beta), ylab="log-likelihood")
lines(xx, pz, lty=2)
lines(xx, pw, lty=3)

