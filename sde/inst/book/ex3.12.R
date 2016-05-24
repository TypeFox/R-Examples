require(sde)

K.est <- function(x) {
  n.obs <- length(x)
  n.obs/(2*(sum(x^2)))
}

LS.est <- function(x) {
  n <- length(x) -1
  k.sum <- sum(x[1:n]*x[2:(n+1)])
  dt <- deltat(x)
  ifelse(k.sum>0, -log(k.sum/sum(x[1:n]^2))/dt, NA)
}

MLE.est <- function(y, lower=0, upper=Inf){
 n <- length(y) - 1
 Dt <- deltat(y)
 Y <- y[2:(n+1)]
 g <- function(theta){
  ss <- sqrt((1-exp(-2*Dt*theta))/(2*theta))
  X <- y[1:n]*exp(-theta*Dt)
  lik <- dnorm(Y,mean=X,sd=ss)  
  -sum(log(lik))
 } 
 tmp <- try(optim(runif(1), g, method="L-BFGS-B", lower=lower, upper=upper)$par) 
 if(class(tmp)=="try-error") tmp <- NA
 tmp
}

# ex3.12.R
set.seed(123); 
d <- expression(-1 * x)
s <- expression(1) 
x0 <- rnorm(1,sd=sqrt(1/2))
sde.sim(X0=x0,drift=d, sigma=s,N=2500,delta=0.1) -> X
 
# Kessler's estimator revisited
f <- list(expression(2*theta*x^2-1))
simple.ef(X, f, lower=0, upper=Inf)

K.est(X)

# Least Squares estimator revisited
f <- list(expression(x*(y-x*exp(-0.1*theta))))
simple.ef(X, f, lower=0, upper=Inf)

LS.est(X)

# MLE estimator revisited, f is based on 
# the score function of the process
f <-list(expression((1 + 2 * (x^2) * theta + 2 * 0.1 * theta + 
   exp(4 * 0.1 * theta) * (1 - 2 * (y^2) * theta) 
   - 4 * exp(3 * 0.1 * theta) * x * y * theta * (-1 + 0.1 * theta)
   - 4 * exp(theta * 0.1) * x * y * theta * (1 + theta * 0.1) 
   + 2 * exp(2 * 0.1 * theta) * (-1 - 0.1 * theta + (x^2) * theta * 
   (-1 + 2 * 0.1 * theta) + (y^2) * theta * (1 + 2 *
   0.1 * theta)))/((2 * (-1 + exp(2 * 0.1 * theta))^2) * theta)))

simple.ef(X, f, lower=0, upper=Inf)

MLE.est(X)


# ex3.12.R (cont)

set.seed(123); 
d <- expression(-1 * x)
s <- expression(1) 
x0 <- rnorm(1,sd=sqrt(1/2))
sde.sim(X0=x0,drift=d, sigma=s,N=1500,delta=0.1) -> X

d <- expression(-theta* x)
s <- expression(1) 
h <- list(expression(x^2))
simple.ef2(X, d, s, h, lower=0, upper=Inf)

K.est(X)
MLE.est(X)
LS.est(X)


# ex3.12.R (cont)

set.seed(123); 
d <- expression(10 - x)
s <- expression(sqrt(x)) 
x0 <- 10
sde.sim(X0=x0,drift=d, sigma=s,N=1500,delta=0.1) -> X

d <- expression(alpha +theta* x)
s <- expression(sqrt(x)) 
h <- list(expression(x),expression(x^2))
simple.ef(X, d, s, h, lower=c(0,-Inf), upper=c(Inf,0))

# explicit estimator for alpha
(sum(X)^2)/(2*(length(X)*sum(X^2)-sum(X)^2))

# explicit estimator for theta
(-length(X)*sum(X))/(2*(length(X)*sum(X^2)-sum(X)^2))


# ex3.12.R (cont)
set.seed(123); 
d <- expression(10 - x)
s <- expression(sqrt(x)) 
x0 <- 10
sde.sim(X0=x0,drift=d, sigma=s,N=1500,delta=0.1) -> X

d <- expression(10- x)
s <- expression(x^gamma) 
h <- list(expression(x^2))
simple.ef2(X, d, s, h, lower=0, upper=1)

