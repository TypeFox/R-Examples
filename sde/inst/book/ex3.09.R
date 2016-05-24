require(sde)
# ex3.09.R

set.seed(123)
d <- expression(-1*x)
s <- expression(2) 
sde.sim(drift=d, sigma=s) -> X

M0 <- function(t, x, theta) -theta[1]*x
M1 <- function(t, x, theta) -theta[1]
M2 <- function(t, x, theta) 0
M3 <- function(t, x, theta) 0
M4 <- function(t, x, theta) 0
M5 <- function(t, x, theta) 0
M6 <- function(t, x, theta) 0
mu <- list(M0, M1, M2, M3, M4, M5, M6)

F <- function(t, x, theta) x/sqrt(theta[2])
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


# ex3.09.R (cont)
xx <- seq(-3,3,length=100)
sapply(xx, function(x) HPloglik(X,c(x,4),mu,F,S)) -> px
sapply(xx, function(x) true.loglik(c(x,4))) -> py
sapply(xx, function(x) EULERloglik(X,c(x,4),B,S)) -> pz

plot(xx,px,type="l",xlab=expression(beta),ylab="log-likelihood") # approx
lines(xx,py, lty=3) # true
lines(xx,pz, lty=2) # Euler


# ex3.09.R (cont)
require(stats4)
HP.negloglik <- function(BETA=3, SIGMA2=2) 
  -HPloglik(X,c(BETA,SIGMA2),mu,F,S)
true.negloglik <- function(BETA=3, SIGMA2=2) 
  -true.loglik(c(BETA,SIGMA2))
euler.negloglik <- function(BETA=3, SIGMA2=2) 
  -EULERloglik(X,c(BETA,SIGMA2),B,S)

mle(true.negloglik,lower=c(0,0),method="L-BFGS-B") -> fit.true
mle(HP.negloglik,lower=c(0,0),method="L-BFGS-B") -> fit.approx
mle(euler.negloglik,lower=c(0,0),method="L-BFGS-B") -> fit.euler

# we look at the estimates
coef(fit.true)
coef(fit.approx)
coef(fit.euler)

# ex3.09.R (cont)
logLik(fit.true)
logLik(fit.approx)
logLik(fit.euler)

# ex3.09.R (cont)
vcov(fit.true)
vcov(fit.approx)
vcov(fit.euler)

# ex3.09.R (cont)
beta <- 1
sigma <- 2
DELTA <- deltat(X)
vbeta <- (exp(2*beta*DELTA)-1)/DELTA^2
cv.bsigma <- sigma^2*(exp(2*beta*DELTA)-1-2*beta*DELTA)/(beta*DELTA^2)
vsigma <- sigma^4 *((exp(2*beta*DELTA)-1)^2+2*beta^2*DELTA^2*(exp(2*beta*DELTA)+1)+4*beta*DELTA*(exp(2*beta*DELTA)-1))/(beta^2*DELTA^2*(exp(2*beta*DELTA)-1))
matrix(c(vbeta, cv.bsigma, cv.bsigma, vsigma),2,2)

# ex3.09.R (cont)
vcov(fit.true)*100 # the sample size

# ex3.09.R (cont)
mle(true.negloglik,lower=c(0,0),fixed=list(SIGMA2=4),method="L-BFGS-B") -> fit.true
mle(HP.negloglik,lower=c(0,0),fixed=list(SIGMA2=4),method="L-BFGS-B") -> fit.approx
mle(euler.negloglik,lower=c(0,0),fixed=list(SIGMA2=4),method="L-BFGS-B") -> fit.euler
coef(fit.true)
coef(fit.approx)
coef(fit.euler)

vcov(fit.true)
vcov(fit.approx)
vcov(fit.euler)

