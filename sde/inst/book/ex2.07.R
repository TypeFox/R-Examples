# ex2.07.R
# Cox-Ingersoll-Ross (CIR-2)
d <- expression( 6-3*x ) 
s <- expression( 2*sqrt(x) ) 
s.x <- expression( 1/sqrt(x) )
require(sde)
set.seed(123)
sde.sim(X0=10, drift=d, sigma=s, sigma.x=s.x, 
    method="milstein") -> X
plot(X,main="Cox-Ingersoll-Ross")
