# ex2.08.R
# Cox-Ingersoll-Ross (CIR-3)
require(sde)
set.seed(123)
d <- expression( (6-3*x^2 - 1)/(2*x) )
s <- expression( 1 )
sde.sim(X0=sqrt(10),drift=d, sigma=s) -> Y
plot(Y^2,main="Cox-Ingersoll-Ross")
