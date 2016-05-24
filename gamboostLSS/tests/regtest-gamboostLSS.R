###
# check gamboostLSS()

require("gamboostLSS")

set.seed(1907)
x1 <- rnorm(1000)
x2 <- rnorm(1000)
x3 <- rnorm(1000)
x4 <- rnorm(1000)
x5 <- rnorm(1000)
x6 <- rnorm(1000)
mu    <- exp(1.5 + 0.3 * x1^2 + 0.5 * x2 - 3 * sin(x3) -1 * x4)
sigma <- exp(-0.2 * x4 + 0.2 * x5 + 0.4 * x6)
y <- numeric(1000)
for( i in 1:1000)
    y[i] <- rnbinom(1, size = sigma[i], mu = mu[i])
dat <- data.frame(x1, x2, x3, x4, x5, x6, y)

model <- gamboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                     control = boost_control(mstop = 200))

coef(model)

par(mfrow = c(3,2))
plot(dat$x3, fitted(model$mu, which = "x3"), main = "mu")
lines(sort(dat$x3), - 3 * sin(dat$x3)[order(dat$x3)], col = "red")
plot(dat$x3, fitted(model$sigma, which = "x3"), main = "sigma")
lines(sort(dat$x3), - 3 * sin(dat$x3)[order(dat$x3)], col = "red")

model[400]
plot(dat$x3, fitted(model$mu, which = "x3"), main = "mu")
lines(sort(dat$x3), - 3 * sin(dat$x3)[order(dat$x3)], col = "red")
plot(dat$x3, fitted(model$sigma, which = "x3"), main = "sigma")
lines(sort(dat$x3), - 3 * sin(dat$x3)[order(dat$x3)], col = "red")

model[600]
plot(dat$x3, fitted(model$mu, which = "x3"), main = "mu")
lines(sort(dat$x3), - 3 * sin(dat$x3)[order(dat$x3)], col = "red")
plot(dat$x3, fitted(model$sigma, which = "x3"), main = "sigma")
lines(sort(dat$x3), - 3 * sin(dat$x3)[order(dat$x3)], col = "red")
