require(sde)

# ex3.16.R
alpha <- 0.5
beta <- 0.2
sigma <- sqrt(0.05)
true <- c(alpha, beta, sigma)
names(true) <- c("alpha", "beta", "sigma")
 
set.seed(123)
x0 <- rsCIR(1,theta=true)
sde.sim(X0=x0,model="CIR",theta=true,N=500000,delta=0.001) -> X
X <- window(X, deltat=0.1)
DELTA = deltat(X)
n <- length(X) 
X <- window(X, start=n*DELTA*0.5)
plot(X)


# ex3.16.R (cont)
u <- function(x, y, theta, DELTA){
  c.mean <- theta[1]/theta[2] + (y-theta[1]/theta[2])*exp(-theta[2]*DELTA)
  c.var <- ((y*theta[3]^2 * 
     (exp(-theta[2]*DELTA)-exp(-2*theta[2]*DELTA))/theta[2] +  
        +theta[1]*theta[3]^2*(1-exp(-2*theta[2]*DELTA))/(2*theta[2]^2)))  
  cbind(x-c.mean,y*(x-c.mean), c.var-(x-c.mean)^2, y*(c.var-(x-c.mean)^2))  
}


# ex3.16.R (cont)
CIR.lik <- function(theta1,theta2,theta3) {
 n <- length(X)
 dt <- deltat(X)
 -sum(dcCIR(x=X[2:n], Dt=dt, x0=X[1:(n-1)], theta=c(theta1,theta2,theta3), 
   log=TRUE))
}

fit <- mle(CIR.lik, start=list(theta1=.1,  theta2=.1,theta3=.3), 
    method="L-BFGS-B",lower=c(0.001,0.001,0.001), upper=c(1,1,1))
# maximum likelihood estimates
coef(fit)

gmm(X,u, guess=as.numeric(coef(fit)), lower=c(0,0,0), upper=c(1,1,1))

true


# ex3.16.R (cont)
u2 <- function(x, y, theta, DELTA){
  c.mean <- theta[1]/theta[2] + (y-theta[1]/theta[2])*exp(-theta[2]*DELTA)
  cbind(x-c.mean,y*(x-c.mean))
}

set.seed(123)
gmm(X, u2, dim=2, lower=c(0,0), upper=c(1,1))
true
