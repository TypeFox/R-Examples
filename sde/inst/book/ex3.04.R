# ex3.04.R

require(stats4)
require(sde)

CIR.lik <- function(theta1,theta2,theta3) {
 n <- length(X)
 dt <- deltat(X)
 -sum(dcCIR(x=X[2:n], t=dt, x0=X[1:(n-1)], theta=c(theta1,theta2,theta3), 
   log=TRUE))
}
 
# non efficient version based on non-central Chi^2 density
dcCIR2 <-function (x, t, x0, theta, log = FALSE) 
{
   c <- 2*theta[2]/((1-exp(-theta[2]*t))*theta[3]^2)
   ncp <- 2*c*x0*exp(-theta[2]*t)
   df <- 4*theta[1]/theta[3]^2
   lik <- (dchisq(2 * x * c, df = df, ncp = ncp, log = TRUE) 
    + log(2*c))
   if(!log)
    lik <- exp(lik) 
  lik
}
CIR.lik2 <- function(theta1,theta2,theta3) {
 n <- length(X)
 dt <- deltat(X)
 -sum(dcCIR2(x=X[2:n], t=dt, x0=X[1:(n-1)], theta=c(theta1,theta2,theta3),
   log=TRUE))
}


set.seed(123)
X <- sde.sim(X0=.1, model="CIR", theta=c(.2, 0.06, 0.15), N=2500, delta=0.1)
 
# ex3.04.R 
system.time(L1 <-CIR.lik(.1,.1,.1))
print(L1, digits=12)
system.time(L2 <- CIR.lik2(.1,.1,.1))
print(L2, digits=12)

 
# ex3.04.R (cont.)
mle(CIR.lik, start=list(theta1=.1,  theta2=.1,theta3=.3), 
  method="L-BFGS-B",lower=c(0.001,0.001,0.001), upper=c(1,1,1)) -> fit
fit
mle(CIR.lik2, start=list(theta1=.1,  theta2=.1,theta3=.3), 
  method="L-BFGS-B",lower=c(0.001,0.001,0.001), upper=c(1,1,1)) -> fit
 