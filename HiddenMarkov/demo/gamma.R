#-----  Test Gamma Distribution  -----
#          estimate rate only

Pi <- matrix(c(0.8, 0.2,
               0.3, 0.7),
             byrow=TRUE, nrow=2)

n <- 1000
x <- dthmm(NULL, Pi, c(0,1), "gamma", pm=list(rate=c(4, 0.5)),
           pn=list(shape=c(rep(3, n/2), rep(5, n/2))))

x <- simulate(x, nsim=n)

#    use above parameter values as initial values
y <- BaumWelch(x)

#   check parameter estimates
print(summary(y))
print(sum(y$delta))
print(y$Pi %*% rep(1, ncol(y$Pi)))


#-----  Test Gamma Distribution  -----
#          estimate shape only

Pi <- matrix(c(0.8, 0.2,
               0.3, 0.7),
             byrow=TRUE, nrow=2)

n <- 1000
x <- dthmm(NULL, Pi, c(0,1), "gamma", pm=list(shape=c(4, 0.1)),
           pn=list(rate=c(rep(0.5, n/2), rep(1, n/2))))

x <- simulate(x, nsim=n)

#    use above parameter values as initial values
y <- BaumWelch(x)

#   check parameter estimates
print(summary(y))
print(sum(y$delta))
print(y$Pi %*% rep(1, ncol(y$Pi)))

