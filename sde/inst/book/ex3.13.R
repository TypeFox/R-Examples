require(sde)
# ex3.13.R
set.seed(123) 
d <- expression(10 - x)
s <- expression(sqrt(x)) 
x0 <- 10
sde.sim(X0=x0,drift=d, sigma=s,N=1500,delta=0.1) -> X

d <- expression(alpha + theta*x)
s <- expression(x^gamma) 
h <- list(expression(x), expression(x^2), expression(x^2))
simple.ef2(X, d, s, h, lower=c(0,-Inf,0), upper=c(Inf,0,1))

# user defined guess
simple.ef2(X, d, s, h, lower=c(0,-Inf,0), upper=c(Inf,0,1), 
   guess=c(1,-.9,.7))

# another guess
simple.ef2(X, d, s, h, lower=c(0,-Inf,0), upper=c(Inf,0,1), 
   guess=c(1,-1.2,.3))


