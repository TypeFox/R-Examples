
library(rockchalk)
N <- 100
dat <- genCorrelatedData(N=N, means=c(100,200), sds=c(20,30), rho=0.4, stde=10)
dat$x3 <- rnorm(100, m=40, s=4)
m1 <- lm(y ~ x1 + x2 + x3, data=dat)
summary(m1)
m1d <- mcDiagnose(m1)

m2 <- lm(y ~ x1 * x2 + x3, data=dat)
summary(m2)
m2d <- mcDiagnose(m2)



m3 <- lm(y ~ log(10+x1) + x3 + poly(x2,2), data=dat)
summary(m3)
m3d <- mcDiagnose(m3)

N <- 100
x1 <- 50 + rnorm(N)
x2 <- log(rgamma(N, 2,1))
x3 <- rpois(N, lambda=17)
z1 <- gl(5, N/5)
dummies <- contrasts(z1)[ as.numeric(z1), ]
dimnames(dummies) <- NULL ## Avoids row name conflict in data.frame below
y3 <- x1  -.5 * x2 + 0.1 * x2^2 + dummies %*% c(0.1,-0.1,-0.2,0.2)+ 5 * rnorm(N)
dat <- data.frame(x1=x1, x2=x2, x3=x3,  z1=z1, y3 = y3)

m3 <- lm(y3 ~ x1 + poly(x2,2)  + log(x1) + z1, dat)
summary(m3)

mcDiagnose(m3)



