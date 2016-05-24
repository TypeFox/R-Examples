#-----  Log Normal Distribution -----

n <- 1000
m <- 2
pm <- list(meanlog=c(3, 4), sdlog=c(0.5, 0.4))

Pi <- matrix(c(0.8, 0.2,
               0.3, 0.7),
             byrow=TRUE, nrow=2)

n <- 1000
x <- dthmm(NULL, Pi, c(1,0), "lnorm", pm=pm)

x <- simulate(x, nsim=n, seed=5)

#    Change initial values for delta
x$delta <- compdelta(Pi)

#    use above parameter values as initial values
y <- BaumWelch(x)

#   check parameter estimates
print(summary(y))
print(sum(y$delta))
print(y$Pi %*% rep(1, ncol(y$Pi)))



if (interactive()) par.default <- par(ask=TRUE)

plot(1:n, x$x, type="l", xlab="Time", ylab="X",
     main=paste("Log Normal ", m, "State HMM"))
for (j in 1:m) points((1:n)[x$y==j], x$x[x$y==j], col=j+1)

w <- seq(0.1, 150, 0.1)
plot(w, dlnorm(w, meanlog=pm$meanlog[1], sdl=pm$sdlog[1]),
     type="l", col=2, ylab="Density",
     main="Log Normal Densities")
points(w, dlnorm(w, meanlog=pm$meanlog[2], sdlog=pm$sdlog[2]),
       type="l", col=3)

error <- residuals(y)

w <- hist(error, main="Log Normal HMM: Pseudo Residuals", xlab="")
z <- seq(-3, 3, 0.01)
points(z, dnorm(z)*n*(w$breaks[2]-w$breaks[1]), col="red", type="l")
box()

qqnorm(error, main="Log Normal HMM: Q-Q Plot of Pseudo Residuals")
abline(a=0, b=1, lty=3)
abline(h=seq(-2, 2, 1), lty=3)
abline(v=seq(-2, 2, 1), lty=3)

if (interactive()) par(par.default)

