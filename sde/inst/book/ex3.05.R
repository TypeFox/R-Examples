# ex3.05.R

require(stats4)
require(sde)

dcEuler <- function(x, t, x0, theta, drift){
 dd <- drift(x0, theta)
 (x-x0)*dd - 0.5*t*dd^2  
}
Euler.lik <- function(theta1,theta2){
  n <- length(X) 
  dt <- deltat(X)
   -sum(dcEuler(X[2:n], dt, X[1:(n-1)], c(theta1,theta2), mydrift))
}
mydrift <- function(x, theta) (theta[1]-theta[2]*x)

set.seed(123)
X <- sde.sim(model="OU", theta =c(5,3,2), N=2500) 
mle(Euler.lik, start=list(theta1=1.5,  theta2=1), 
    method="L-BFGS-B", lower=c(0,0)) -> fit
  
summary(fit)
sqrt(mean((X[2:length(X)] - X[1:(length(X)-1)])^2)/deltat(X))

X <- sde.sim(model="OU", theta =c(5,3,2), delta=0.1, N=2500) 
mle(Euler.lik, start=list(theta1=1.5,  theta2=1), 
    method="L-BFGS-B", lower=c(0,0)) -> fit
  
summary(fit)
sqrt(mean((X[2:length(X)] - X[1:(length(X)-1)])^2)/deltat(X))

# ex3.05.R (cont.)
mle(OU.lik, start=list(theta1=1.5, theta2=1, theta3=1), 
  method="L-BFGS-B", lower=c(0,0,0)) -> fit2
summary(fit2)

