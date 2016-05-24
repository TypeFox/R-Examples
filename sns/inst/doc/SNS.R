### R code from vignette source 'SNS.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: SNS.Rnw:109-110
###################################################
options(prompt = "R> ", continue = "+  ", width = 80, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: SNS.Rnw:313-316
###################################################
library("sns")
library("mvtnorm")
my.seed <- 0


###################################################
### code chunk number 3: SNS.Rnw:319-326
###################################################
logdensity.mvg <- function(x, mu, isigsq) {
  f <- dmvnorm(x = as.numeric(x),
    mean = mu, sigma = solve(isigsq), log = TRUE)
  g <- - isigsq %*% (x - mu)
  h <- -isigsq
  return (list(f = f, g = g, h = h))
}


###################################################
### code chunk number 4: SNS.Rnw:329-338
###################################################
set.seed(my.seed)
K <- 3
mu <- runif(K, min = -0.5, max = +0.5)
isigsq <- matrix(runif(K*K, min = 0.1, max = 0.2), ncol = K)
isigsq <- 0.5*(isigsq + t(isigsq))
diag(isigsq) <- rep(0.5, K)
x.init <- rep(0.0, K)
x.smp <- sns.run(x.init, logdensity.mvg, niter = 500,
  mh.diag = TRUE, mu = mu, isigsq = isigsq)


###################################################
### code chunk number 5: SNS.Rnw:341-342
###################################################
summary(x.smp)


###################################################
### code chunk number 6: SNS.Rnw:352-357
###################################################
library("RegressionFactory")
loglike.poisson <- function(beta, X, y) {
  regfac.expand.1par(beta, X = X, y = y,
    fbase1 = fbase1.poisson.log)
}


###################################################
### code chunk number 7: SNS.Rnw:360-366
###################################################
set.seed(my.seed)
K <- 5
N <- 1000
X <- matrix(runif(N * K, -0.5, +0.5), ncol = K)
beta <- runif(K, -0.5, +0.5)
y <- rpois(N, exp(X %*% beta))


###################################################
### code chunk number 8: SNS.Rnw:369-372
###################################################
beta.init <- rep(0.0, K)
beta.glm <- glm(y ~ X - 1, family = "poisson",
  start = beta.init)$coefficients


###################################################
### code chunk number 9: SNS.Rnw:375-379
###################################################
beta.sns <- sns.run(beta.init, fghEval = loglike.poisson,
  niter = 20, nnr = 20, X = X, y = y)
beta.nr <- beta.sns[20, ]
cbind(beta.glm, beta.nr)


###################################################
### code chunk number 10: SNS.Rnw:382-384
###################################################
beta.smp <- sns.run(beta.init, loglike.poisson
  , niter = 200, nnr = 20, mh.diag = TRUE, X = X, y = y)


###################################################
### code chunk number 11: lp_plot
###################################################
plot(beta.smp, select = 1)


###################################################
### code chunk number 12: fig1
###################################################
plot(beta.smp, select = 1)


###################################################
### code chunk number 13: SNS.Rnw:398-399
###################################################
summary(beta.smp)


###################################################
### code chunk number 14: SNS.Rnw:404-407
###################################################
beta.smp <- sns.run(beta.init, loglike.poisson,
  niter = 1000, nnr = 20, mh.diag = TRUE, X = X, y = y)
predmean.poisson <- function(beta, Xnew) exp(Xnew %*% beta)


###################################################
### code chunk number 15: SNS.Rnw:410-412
###################################################
ymean.new <- predict(beta.smp, predmean.poisson,
  nburnin = 100, Xnew = X)


###################################################
### code chunk number 16: SNS.Rnw:417-421
###################################################
predsmp.poisson <- function(beta, Xnew)
  rpois(nrow(Xnew), exp(Xnew %*% beta))
ysmp.new <- predict(beta.smp, predsmp.poisson
  , nburnin = 100, Xnew = X)


###################################################
### code chunk number 17: SNS.Rnw:424-425
###################################################
summary(ymean.new)


###################################################
### code chunk number 18: SNS.Rnw:431-432
###################################################
summary(ysmp.new)


###################################################
### code chunk number 19: SNS.Rnw:444-453
###################################################
set.seed(my.seed)
K <- 100
X <- matrix(runif(N * K, -0.5, +0.5), ncol = K)
beta <- runif(K, -0.5, +0.5)
y <- rpois(N, exp(X %*% beta))
beta.init <- glm(y ~ X - 1, family = "poisson")$coefficients
beta.smp <- sns.run(beta.init, loglike.poisson,
  niter = 100, nnr = 10, mh.diag = TRUE, X = X, y = y)
summary(beta.smp)


###################################################
### code chunk number 20: SNS.Rnw:461-465
###################################################
beta.smp.part <- sns.run(beta.init, loglike.poisson,
  niter = 100, nnr = 10, mh.diag = TRUE,
  part = sns.make.part(K, 10), X = X, y = y)
summary(beta.smp.part)


###################################################
### code chunk number 21: ssp_plot
###################################################
par(mfrow = c(1,2))
plot(beta.smp, select = 1)
plot(beta.smp.part, select = 1)


###################################################
### code chunk number 22: fig1
###################################################
par(mfrow = c(1,2))
plot(beta.smp, select = 1)
plot(beta.smp.part, select = 1)


###################################################
### code chunk number 23: SNS.Rnw:494-506
###################################################
loglike.linreg.het <- function(coeff, X, Z, y) {
  K1 <- ncol(X)
  K2 <- ncol(Z)
  beta <- coeff[1:K1]
  gamma <- coeff[K1 + 1:K2]
  
  mu <- X %*% beta
  sigma <- sqrt(exp(Z %*% gamma))
  f <- sum(dnorm(y, mu, sigma, log = TRUE))
  
  return (f)
}


###################################################
### code chunk number 24: SNS.Rnw:509-520
###################################################
set.seed(my.seed)
K1 <- 5
K2 <- 5
N <- 1000
X <- matrix(runif(N * K1, -0.5, +0.5), ncol = K1)
Z <- matrix(runif(N * K2, -0.5, +0.5), ncol = K2)
beta <- runif(K1, -0.5, +0.5)
gamma <- runif(K1, -0.5, +0.5)
mu <- X %*% beta
var <- exp(Z %*% gamma)
y <- rnorm(N, X %*% beta, sd = sqrt(var))


###################################################
### code chunk number 25: SNS.Rnw:523-528
###################################################
coeff.init <- rep(0.0, K1 + K2)
check.logdensity <- sns.check.logdensity(coeff.init, loglike.linreg.het
  , X = X, Z = Z, y = y, dx = 1, nevals = 10
  , blocks = list(1:(K1+K2), 1:K1, K1 + 1:K2))
check.logdensity


###################################################
### code chunk number 26: SNS.Rnw:552-566
###################################################
loglike.linreg.het.beta <- function(beta, gamma, X, Z, y) {
  K1 <- length(beta)
  ret <- regfac.expand.2par(c(beta, gamma), X, Z, y
    , fbase2 = fbase2.gaussian.identity.log)
  return (list(f = ret$f, g = ret$g[1:K1], h = ret$h[1:K1, 1:K1]))
}
loglike.linreg.het.gamma <- function(gamma, beta, X, Z, y) {
  K1 <- length(beta)
  K2 <- length(gamma)
  ret <- regfac.expand.2par(c(beta, gamma), X, Z, y
    , fbase2 = fbase2.gaussian.identity.log)
  return (list(f = ret$f, g = ret$g[K1 + 1:K2]
           , h = ret$h[K1 + 1:K2, K1 + 1:K2]))
}


###################################################
### code chunk number 27: SNS.Rnw:569-586
###################################################
nsmp <- 100
beta.iter <- rep(0.0, K1)
gamma.iter <- rep(0.0, K2)
beta.smp <- array(NA, dim = c(nsmp, K1))
gamma.smp <- array(NA, dim = c(nsmp, K1))
for (n in 1:nsmp) {
  beta.iter <- sns(beta.iter, loglike.linreg.het.beta
    , gamma = gamma.iter, X = X, Z = Z, y = y, rnd = nsmp>10)
  gamma.iter <- sns(gamma.iter, loglike.linreg.het.gamma
    , beta = beta.iter, X = X, Z = Z, y = y, rnd = nsmp>10)
  beta.smp[n, ] <- beta.iter
  gamma.smp [n, ] <- gamma.iter
}
beta.est <- colMeans(beta.smp[(nsmp/2+1):nsmp, ])
gamma.est <- colMeans(gamma.smp[(nsmp/2+1):nsmp, ])
cbind(beta, beta.est)
cbind(gamma, gamma.est)


###################################################
### code chunk number 28: SNS.Rnw:589-608
###################################################
library(MfUSampler)
loglike.linreg.het.gamma.fonly <- function(gamma, beta, X, Z, y) {
  return (regfac.expand.2par(c(beta, gamma), X, Z, y
    , fbase2 = fbase2.gaussian.identity.log, fgh = 0))
}
beta.iter <- rep(0.0, K1)
gamma.iter <- rep(0.0, K2)
for (n in 1:nsmp) {
  beta.iter <- sns(beta.iter, loglike.linreg.het.beta
    , gamma = gamma.iter, X = X, Z = Z, y = y, rnd = nsmp>10)
  gamma.iter <- MfU.Sample(gamma.iter, loglike.linreg.het.gamma.fonly
    , beta = beta.iter, X = X, Z = Z, y = y)
  beta.smp[n, ] <- beta.iter
  gamma.smp [n, ] <- gamma.iter
}
beta.est.mix <- colMeans(beta.smp[(nsmp/2+1):nsmp, ])
gamma.est.mix <- colMeans(gamma.smp[(nsmp/2+1):nsmp, ])
cbind(beta, beta.est.mix)
cbind(gamma, gamma.est.mix)


###################################################
### code chunk number 29: SNS.Rnw:612-628
###################################################
loglike.linreg.het.gamma.numaug <- 
  sns.fghEval.numaug(loglike.linreg.het.gamma.fonly, numderiv = 2)
beta.iter <- rep(0.0, K1)
gamma.iter <- rep(0.0, K2)
for (n in 1:nsmp) {
  beta.iter <- sns(beta.iter, loglike.linreg.het.beta
    , gamma = gamma.iter, X = X, Z = Z, y = y, rnd = nsmp>10)
  gamma.iter <- sns(gamma.iter, loglike.linreg.het.gamma
    , beta = beta.iter, X = X, Z = Z, y = y, rnd = nsmp>10)
  beta.smp[n, ] <- beta.iter
  gamma.smp [n, ] <- gamma.iter
}
beta.est.num <- colMeans(beta.smp[(nsmp/2+1):nsmp, ])
gamma.est.num <- colMeans(gamma.smp[(nsmp/2+1):nsmp, ])
cbind(beta, beta.est.num)
cbind(gamma, gamma.est.num)


###################################################
### code chunk number 30: SNS.Rnw:697-698
###################################################
sessionInfo()


