# ex3.02.R
OU.lik <- function(theta1, theta2, theta3){
  n <- length(X)
  dt <- deltat(X)
  -sum(dcOU(X[2:n], dt, X[1:(n-1)], c(theta1,theta2,theta3), log=TRUE))
 }

require(stats4)
require(sde)
set.seed(123)
X <- sde.sim(model="OU", theta=c(0,3,2), N=1000, delta=1)
mle(OU.lik, start=list(theta2=1.5, theta3=1), fixed=list(theta1=0), 
  method="L-BFGS-B", lower=c(0,0)) -> fit
summary(fit)

# ex3.02.R (cont.)
n <- length(X) 
tmp.sum <- sum(X[1:(n-1)]*X[2:n])
dt <- deltat(X)
theta2.hat <- ifelse(tmp.sum>0, -log(tmp.sum/sum(X[1:(n-1)]^2))/dt ,NA)
theta2.hat
theta3sq.hat <- 2*theta2.hat/((n-1)*(1-exp(-2*dt*theta2.hat))) * 
     sum((X[2:n]-X[1:(n-1)]*exp(-dt*theta2.hat))^2)
sqrt(theta3sq.hat)


# ex3.02.R (cont.)
set.seed(123)
X <- sde.sim(model="OU", theta=c(0,3,1), N=1000, delta=1)
mle(OU.lik, start=list(theta2=1.5), 
      method="L-BFGS-B", lower=c(0), fixed=c(theta1=0, theta3=1)) -> fit
summary(fit)
