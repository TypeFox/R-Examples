require(sde)

LS.est <- function(x) {
  n <- length(x) -1
  k.sum <- sum(x[1:n]*x[2:(n+1)])
  dt <- deltat(x)
  ifelse(k.sum>0, -log(k.sum/sum(x[1:n]^2))/dt, NA)
}

# ex3.14.R
set.seed(123)
d <- expression(-1 * x)
s <- expression(1) 
x0 <- rnorm(1,sd=sqrt(1/2))
sde.sim(X0=x0,drift=d, sigma=s,N=1000,delta=0.1) -> X
 
d <- expression(-theta * x)
  
linear.mart.ef(X, d, s, a1=expression(-x), lower=0, upper=Inf,
  c.mean=expression(x*exp(-theta*0.1)), 
  c.var=expression((1-exp(-2*theta*0.1))/(2*theta)))


# the linear mart. e.f. coincides with the
# least square estimator
LS.est(X)

