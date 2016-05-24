# ex2.09.R
# geometric Brownian Motion
require(sde)
set.seed(123)
d <- expression( x ) 
s <- expression( 0.5*x ) 
sde.sim(X0=10,drift=d, sigma=s) -> X
plot(X,main="geometric Brownian Motion")
