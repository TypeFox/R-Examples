# ex2.12.R
set.seed(123)
X <- sde.sim(drift=expression(-3*x), method="ozaki") 
set.seed(123)
Y <- sde.sim(drift=expression(-3*x))
plot(X)
lines(as.numeric(time(Y)), Y, col="red")
