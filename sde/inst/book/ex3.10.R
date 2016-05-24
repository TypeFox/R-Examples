require(sde)
# ex3.10.R

F <- function(t, x, theta) 2*sqrt(x)/theta[3]
S <- function(t, x, theta) theta[3]*sqrt(x)

M0 <- function(t, x, theta) 
 (4*theta[1]*theta[2]-theta[3]^2)/(2*x*theta[3]^2) - 0.5*theta[1]*x
M1 <- function(t, x, theta) 
 -0.5*theta[1]-(4*theta[1]*theta[2]-theta[3]^2)/(2*x^2*theta[3]^2)
M2 <- function(t, x, theta) 
 (4*theta[1]*theta[2]-theta[3]^2)/(x^3*theta[3]^2)
M3 <- function(t, x, theta) 
 -3*(4*theta[1]*theta[2]-theta[3]^2)/(x^4*theta[3]^2)
M4 <- function(t, x, theta) 
 12*(4*theta[1]*theta[2]-theta[3]^2)/(x^5*theta[3]^2)
M5 <- function(t, x, theta) 
 -60*(4*theta[1]*theta[2]-theta[3]^2)/(x^6*theta[3]^2)
M6 <- function(t, x, theta) 
 360*(4*theta[1]*theta[2]-theta[3]^2)/(x^7*theta[3]^2)
mu1 <- list(M0, M1, M2, M3, M4, M5, M6)

# ex3.10.R (cont)
# we now ask R to calculate derivatives
m0 <- expression((4*theta1*theta2-theta3^2)/(2*x*theta3^2) - 0.5*theta1*x)

params <- all.vars(m0)
params <- params[-which(params=="x")]
np <- length(params)

# we construct derivatives by iteration
for(i in 1:6){
 esp <- get(sprintf("m%d",i-1))
 assign(sprintf("m%d",i), D(esp, "x"))
}

mu2 <- vector(7, mode="list")
# `mu2' must be a list of functions in (t,x,theta)
mu2[[1]] <- function(t,x,theta) eval(m0)   
mu2[[2]] <- function(t,x,theta) eval(m1)   
mu2[[3]] <- function(t,x,theta) eval(m2)   
mu2[[4]] <- function(t,x,theta) eval(m3)   
mu2[[5]] <- function(t,x,theta) eval(m4)   
mu2[[6]] <- function(t,x,theta) eval(m5)   
mu2[[7]] <- function(t,x,theta) eval(m6)   


# we need to remap theta[1:3] into (theta1,theta2,theta3)
# hence we write a wrapper function that calls HPloglik
# in the correct way
HPloglik2 <- function(X, theta, mu, F, S){
 sapply(1:np, function(x) assign(params[x], theta[x], .GlobalEnv))
 HPloglik(X, theta, mu2, F, S)
}

# ex3.10.R (cont)
set.seed(123)
X <- sde.sim(X0=1, model="CIR", theta=c(.2, .4, .15),delta=1e-3,N=100)

xx <- seq(0,4,length=100)
a <- system.time(sapply(xx, function(x) HPloglik(X, c(.4,x/.4,.15), mu1, F, S))-> px1)
b <- system.time(sapply(xx, function(x) HPloglik2(X, c(.4,x/.4,.15), mu2, F, S))-> px2)

# should be zero
sum(abs(px1-px2))

# should be very high
b/a

# ex3.10.R (cont)
# true CIR log-likelihood
CIR.lik <- function(alpha,kappa,sigma) {
 n <- length(X)
 dt <- deltat(X)
 sum(dcCIR(x=X[2:n], Dt=dt, x0=X[1:(n-1)], theta=c(alpha,kappa,sigma), 
   log=TRUE))
}

CIR.negloglik <- function(ALPHA=1, KAPPA=1, SIGMA=1) 
  -CIR.lik(ALPHA,KAPPA,SIGMA)

HP.negloglik <- function(THETA1=1, THETA2=1, THETA3=1) 
  -HPloglik(X,c(THETA1,THETA2,THETA3),mu1,F,S)

mle(HP.negloglik,lower=c(0.01,0.01,0.01),fixed=list(THETA1=.4,THETA3=.15),
 method="L-BFGS-B") -> fit.approx
THETA <- coef(fit.approx)
THETA
# alpha = theta1*theta2
cat(sprintf("alpha=%f\n",THETA[1]*THETA[2]))

mle(CIR.negloglik,lower=c(0.01,0.01,0.01),fixed=list(KAPPA=.4,SIGMA=.15),
 method="L-BFGS-B") -> fit.CIR
PAR <- coef(fit.CIR)
PAR

xx <- seq(0.0858,0.086,length=100)
sapply(xx, function(x) HPloglik(X, c(.4,x/.4,.15), mu1, F, S))-> px1
sapply(xx, function(x) CIR.lik(x,.4,.15))-> px3

par(mar=c(5,5,1,1))
plot(xx,px1,type="l",xlab=expression(alpha),ylab="log-likelihood")
lines(xx,px3,lty=2)
abline(v=PAR[1],lty=3)
abline(v=prod(THETA[1:2]),lty=3)

# ex3.10.R (cont)
theta <- c(.2, .4, .15)
sapply(1:3, function(x) assign(params[x], theta[x], .GlobalEnv))
mu1[[4]]
mu2[[4]]
mu1[[4]](0,1,theta)
mu2[[4]](0,1,theta)

