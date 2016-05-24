if (interactive()) par.default <- par(ask=TRUE)

#-----  Test Logistic Distribution  -----

Pi <- matrix(c(0.8, 0.2,
               0.3, 0.7),
             byrow=TRUE, nrow=2)

pm <- list(location=c(8, -2), scale=c(1, 0.5))

x <- dthmm(NULL, Pi, c(0,1), "logis", pm=pm)

x <- simulate(x, nsim=1000)

hist(x$x, main="", xlab=expression(x))
box()

#    use above parameter values as initial values
y <- BaumWelch(x)

#   check parameter estimates
print(summary(y))
print(sum(y$delta))
print(y$Pi %*% rep(1, ncol(y$Pi)))


#---------------------------------------------
#    Fixed Scale Parameter

Pi <- matrix(c(0.8, 0.2,
               0.3, 0.7),
             byrow=TRUE, nrow=2)

n <- 1000
pm <- list(location=c(8, -2))
pn <- list(scale=rep(1, n))

x <- dthmm(NULL, Pi, c(0,1), "logis", pm=pm, pn=pn)

x <- simulate(x, nsim=n)

hist(x$x, main="", xlab=expression(x))
box()

#    use above parameter values as initial values
y <- BaumWelch(x)

#   check parameter estimates
print(summary(y))
print(sum(y$delta))
print(y$Pi %*% rep(1, ncol(y$Pi)))

z <- residuals(y)
qqnorm(z, main="Logistic HMM: Q-Q Plot of Pseudo Residuals")
abline(a=0, b=1, lty=3)
abline(h=seq(-2, 2, 1), lty=3)
abline(v=seq(-2, 2, 1), lty=3)


#---------------------------------------------
#    Fixed Location Parameter

Pi <- matrix(c(0.8, 0.2,
               0.3, 0.7),
             byrow=TRUE, nrow=2)

n <- 1000
pm <- list(scale=c(0.5, 2))
pn <- list(location=c(rep(5, n/2), rep(1, n/2)))

x <- dthmm(NULL, Pi, c(0,1), "logis", pm=pm, pn=pn)

x <- simulate(x, nsim=n)

hist(x$x, main="", xlab=expression(x))
box()

#    use above parameter values as initial values
y <- BaumWelch(x)

#   check parameter estimates
print(summary(y))
print(sum(y$delta))
print(y$Pi %*% rep(1, ncol(y$Pi)))

z <- residuals(y)
qqnorm(z, main="Logistic HMM: Q-Q Plot of Pseudo Residuals")
abline(a=0, b=1, lty=3)
abline(h=seq(-2, 2, 1), lty=3)
abline(v=seq(-2, 2, 1), lty=3)

if (interactive()) par(par.default)
