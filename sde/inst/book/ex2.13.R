# ex2.13.R
bX <- expression((5 - 11 * x + 6 * x^2 - x^3))
x0 <- 5
DT <- 0.1
par(mfrow=c(2,3))
set.seed(123)
X <- sde.sim(drift=bX, delta=DT,X0=x0)
plot(X,main="Euler")
set.seed(123)
Y <- sde.sim(drift=bX, method="ozaki",delta=DT,X0=x0) 
plot(Y,main="Ozaki")
set.seed(123)
Z <- sde.sim(drift=bX, method="shoji",delta=DT,X0=x0) 
plot(Z,main="Shoji-Ozaki")
 
DT <- 0.25
set.seed(123)
X <- sde.sim(drift=bX, delta=DT,X0=x0)
plot(X, main="Euler")
set.seed(123)
Y <- sde.sim(drift=bX, method="ozaki",delta=DT,X0=x0) 
plot(Y,main="Ozaki")
set.seed(123)
Z <- sde.sim(drift=bX, method="shoji",delta=DT,X0=x0) 
plot(Z,main="Shoji-Ozaki")
