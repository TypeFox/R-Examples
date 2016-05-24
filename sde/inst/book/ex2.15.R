# ex2.15.R
set.seed(123)
d <- expression(sin(x))
d.x <- expression(cos(x)) 
A <- function(x) 1-cos(x)
sde.sim(method="EA", delta=1/20, X0=0, N=500, drift=d, drift.x = d.x, A=A) -> X
plot(X, main="Periodic drift")

