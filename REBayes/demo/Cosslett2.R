# Cosslett (1983) nonparametric binary response estimator

#eg2: F  <-  1/2 \delta_{-1} + 1/2 \delta_{1} Two components
F  <-  function(x) 0 * (x <= -1) + 1/2 * (x >= -1) + 1/2 * (x >= 1)
n  <-  5000
v  <-  rnorm(n)
eta  <-  sample(c(-1,1),n,replace=TRUE,prob = c(1/2,1/2))
y  <-  (eta < v)-0
fit  <-  Cosslett(x = v, y = y)
par(mfrow = c(1,2))
plot(fit$x,fit$y/sum(fit$y),type = "l", xlab = "x", ylab = "f(x)", main = "density")
plot(fit$x,cumsum(fit$y)/(sum(fit$y)),type = "l", xlab = "x", 
     ylab = "f(x)", main = "distribution F")
lines(fit$x,F(fit$x),col = 3)

