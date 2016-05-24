require(sde)

# ex3.15.R
k <- NULL
theta <- NULL
k1 <- NULL
theta1 <- NULL
sigma.hat <- NULL
sigma <- 0.15
pars <- c(0.5, 0.2, sigma)
n.sim <- 1000
Delta <- 0.01
set.seed(123)
x0 <- rsCIR(n.sim, pars)
 
T <- 100

for(i in 1:n.sim){
  X <- sde.sim(X0=x0[i], model="CIR", theta=pars, N=T/Delta, delta=Delta)
  
  n <- length(X)

# CIR
 I3 <- Delta * sum(1/X[1:(n-1)] + 1/X[2:n])/2
 I1 <- log(X[n]/X[1]) + 0.5*sigma^2 * I3
 I2 <- X[n] - X[1]
 I4 <- Delta * sum(X[1:(n-1)] + X[2:n])/2
 I5 <- n*Delta

 k <- c(k, (I1*I5-I2*I3)/(I3*I4-I5^2))
 theta <- c(theta, (I1*I4-I2*I5)/(I1*I5-I2*I3))
 sigma.est <- sqrt(sum((X[2:n]-X[1:(n-1)])^2/X[1:(n-1)])/(n*Delta))
 sigma.hat <- c(sigma.hat, sigma.est)

 I1 <- log(X[n]/X[1]) + 0.5*sigma.est^2 * I3
 k1 <- c(k1, (I1*I5-I2*I3)/(I3*I4-I5^2))
 theta1 <- c(theta1, (I1*I4-I2*I5)/(I1*I5-I2*I3))

}
cat(sprintf("kappa=%f, theta=%f, sigma=%f : kappa1=%f, theta1=%f\n", mean(k*theta), 
 mean(k), mean(sigma.hat), mean(k1*theta1), mean(k1)))

cat(sprintf("SD: kappa=%f, theta=%f, sigma=%f : kappa1=%f, theta1=%f\n", sd(k*theta), 
 sd(k), sd(sigma.hat), sd(k1*theta1), sd(k1)))

cat(sprintf("kappa=%f, theta=%f, sigma=%f\n", pars[1], pars[2], pars[3]))

