# Cosslett (1983) nonparametric binary response estimator

#eg1: F  <-  N(0,1)
set.seed(12)
n  <-  5000
v  <-  rnorm(n)
eta <- rnorm(n)
y  <-  (eta < v)-0
fit  <-  Cosslett(x = v, y = y)
par(mfrow = c(1,2))
plot(fit$x,fit$y/sum(fit$y),type = "l",xlab = "x", ylab = "f(x)", main = "density")
plot(fit$x,cumsum(fit$y)/(sum(fit$y)),type = "l", xlab = "x", 
     ylab = "F(x)",main = "distribution F")
lines(fit$x,pnorm(fit$x),col = 3)
