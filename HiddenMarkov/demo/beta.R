if (interactive()) par.default <- par(ask=TRUE)

#-----  Test Beta Distribution  -----
#    "shape1" Markov dependent
#    "shape2" time dependent

Pi <- matrix(c(0.8, 0.2,
               0.3, 0.7),
             byrow=TRUE, nrow=2)

x <- seq(0.01, 0.99, 0.01)
plot(x, dbeta(x, shape1=0.5, shape2=2), type="l",
     col="blue", ylab="Density")
points(x, dbeta(x, shape1=8, shape2=2), type="l", col="red")

n <- 1000
x <- dthmm(NULL, Pi, c(0,1), "beta", pm=list(shape1=c(0.5, 6)),
           pn=list(shape2=rep(2, n)))

x <- simulate(x, nsim=n)

#    use above parameter values as initial values
y <- BaumWelch(x)

#   check parameter estimates
print(summary(y))
print(sum(y$delta))
print(y$Pi %*% rep(1, ncol(y$Pi)))


#-----  Test Beta Distribution  -----
#    "shape2" Markov dependent
#    "shape1" time dependent

Pi <- matrix(c(0.8, 0.2,
               0.3, 0.7),
             byrow=TRUE, nrow=2)

x <- seq(0.01, 0.99, 0.01)
plot(x, dbeta(x, shape1=2, shape2=6), type="l",
     col="blue", ylab="Density")
points(x, dbeta(x, shape1=2, shape2=0.5), type="l", col="red")

n <- 1000
x <- dthmm(NULL, Pi, c(0,1), "beta", pm=list(shape2=c(0.5, 6)),
           pn=list(shape1=rep(2, n)))

x <- simulate(x, nsim=n)

#    use above parameter values as initial values
y <- BaumWelch(x)

#   check parameter estimates
print(summary(y))
print(sum(y$delta))
print(y$Pi %*% rep(1, ncol(y$Pi)))

if (interactive()) par(par.default)

