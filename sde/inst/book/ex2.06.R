# ex2.06.R
# Cox-Ingersoll-Ross (CIR-1)
require(sde)
set.seed(123)
d <- expression( 6-3*x ) 
s <- expression( 2*sqrt(x) ) 
sde.sim(X0=10,drift=d, sigma=s) -> X
plot(X,main="Cox-Ingersoll-Ross")
