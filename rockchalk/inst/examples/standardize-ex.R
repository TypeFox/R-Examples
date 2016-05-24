
library(rockchalk)
N <- 100
dat <- genCorrelatedData(N = N, means = c(100,200), sds = c(20,30), rho = 0.4, stde = 10)
dat$x3 <- rnorm(100, m = 40, s = 4)

m1 <- lm(y ~ x1 + x2 + x3, data = dat)
summary(m1)

m1s <- standardize(m1)
summary(m1s)



m2 <- lm(y ~ x1 * x2 + x3, data = dat)
summary(m2)

m2s <- standardize(m2)
summary(m2s)

m2c <- meanCenter(m2)
summary(m2c)

