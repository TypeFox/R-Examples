# ex3.03.R
dcBS <- function(x, t, x0, theta, log = TRUE){
   ml <- log(x0) + (theta[1]-theta[2]^2/2)*t
   sl <- sqrt(t)*theta[2]
   lik <- dlnorm(x, meanlog = ml, sdlog = sl, log=TRUE)
  if(!log)
   lik <- exp(lik)
  lik
}
BS.lik <- function(theta1,theta2) {
 n <- length(X)
 dt <- deltat(X)
 -sum(dcBS(x=X[2:n], t=dt, x0=X[1:(n-1)], theta=c(theta1,theta2),
  log=TRUE))
}

require(stats4)
require(sde)


set.seed(123)
X <- sde.sim(model="BS", theta=c(.5,.2), delta=0.01)
mle(BS.lik, start=list(theta1=1,  theta2=1), 
    method="L-BFGS-B", lower=c(0.01,0.01)) -> fit
coef(fit)
length(X)*deltat(X)

set.seed(123)
X <- sde.sim(model="BS", theta=c(.5,.2), N=1000, delta=0.01)
mle(BS.lik, start=list(theta1=1,  theta2=1), 
    method="L-BFGS-B", lower=c(0.01,0.01)) -> fit
coef(fit)
length(X)*deltat(X)

set.seed(123)
X <- sde.sim(model="BS", theta=c(.5,.2), N=5000, delta=0.01)
mle(BS.lik, start=list(theta1=1,  theta2=1), 
    method="L-BFGS-B", lower=c(0.01,0.01)) -> fit
coef(fit)
length(X)*deltat(X)
