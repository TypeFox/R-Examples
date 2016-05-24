# ex3.01.R
OU.lik <- function(theta1, theta2, theta3){
  n <- length(X)
  dt <- deltat(X)
  -sum(dcOU(X[2:n], dt, X[1:(n-1)], c(theta1,theta2,theta3), log=TRUE))
 }

require(stats4)
require(sde)
set.seed(123)
X <- sde.sim(model="OU", theta=c(3,1,2), N=1000, delta=1)
mle(OU.lik, start=list(theta1=1, theta2=0.5, theta3=1), 
      method="L-BFGS-B", lower=c(-Inf,0,0)) -> fit
summary(fit)

# ex3.01.R (cont.)
prof <- profile(fit)
par(mfrow=c(1,3))
plot(prof)
par(mfrow=c(1,1))
vcov(fit)
confint(fit)
 
# ex3.01.R (cont.)
set.seed(123)
X <- sde.sim(model="OU", theta=c(3,1,2), N=1000, delta=1e-3)
mle(OU.lik, start=list(theta1=1, theta2=0.5, theta3=1), 
      method="L-BFGS-B", lower=c(-Inf,0,0)) -> fit2
summary(fit2)
