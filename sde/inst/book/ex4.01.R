# ex4.01.R
set.seed(123)

dri <- expression(-(x-10))
dif <- expression(2*sqrt(x)) 
sde.sim(X0=10,drift=dri, sigma=dif,N=1000,delta=0.1) -> X

b <- function(x,theta) -theta[1]*(x-theta[2])
b.x <- function(x,theta)  -theta[1]+0*x

s <- function(x,theta) theta[3]*sqrt(x)
s.x <- function(x,theta) theta[3]/(2*sqrt(x))
s.xx <- function(x,theta) -theta[3]/(4*x^1.5)

# we let sdeAIC calculate the estimates and the AIC statistics
sdeAIC(X, NULL, b, s, b.x, s.x, s.xx, guess=c(1,1,1),
            lower=rep(1e-3,3), method="L-BFGS-B")


