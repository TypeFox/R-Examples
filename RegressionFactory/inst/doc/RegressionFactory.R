### R code from vignette source 'RegressionFactory.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: RegressionFactory.Rnw:101-102
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: RegressionFactory.Rnw:232-245 (eval = FALSE)
###################################################
## regfac.expand.1par <- function(beta, X, y, fbase1, fgh = 2, ...) {
##   # obtain base distribution derivatives
##   ret <- fbase1(X %*% beta, y, fgh, ...)
##   # expand base derivatives
##   f <- sum(ret$f)
##   if (fgh == 0) return (f)
##   g <- t(X) %*% ret$g
##   if (fgh == 1) return (list(f = f, g = g))
##   xtw <- 0*X
##   for (k in 1:ncol(X)) xtw[, k] <- X[, k] * ret$h
##   h <- t(xtw) %*% X
##   return (list(f = f, g = g, h = h))
## }


###################################################
### code chunk number 3: RegressionFactory.Rnw:256-295 (eval = FALSE)
###################################################
## regfac.expand.2par <- function(coeff, X
##   , Z=matrix(1.0, nrow = nrow(X), ncol = 1)
##   , y, fbase2, fgh = 2, block.diag = FALSE
##   , ...) {
##   # extracting coefficients of X and Z
##   K1 <- ncol(X); K2 <- ncol(Z)
##   beta <- coeff[1:K1]
##   gamma <- coeff[K1 + 1:K2]
##   
##   # obtain base distribution derivatives
##   ret <- fbase2(X %*% beta, Z %*% gamma, y, fgh, ...)
## 
##   # expand base derivatives
##   # function
##   f <- sum(ret$f)
##   if (fgh == 0) return (f)
##   # gradient
##   g <- c(t(X) %*% ret$g[, 1], t(Z) %*% ret$g[, 2])
##   if (fgh == 1) return (list(f = f, g = g))
##   # Hessian
##   h <- array(0, dim=c(K1+K2, K1+K2))
##   # XX block
##   xtw <- 0 * X
##   for (k in 1:K1) xtw[, k] <- X[, k] * ret$h[, 1]
##   h[1:K1, 1:K1] <- t(xtw) %*% X
##   # ZZ block
##   ztw <- 0 * Z
##   for (k in 1:K2) ztw[, k] <- Z[, k] * ret$h[, 2]
##   h[K1 + 1:K2, K1 + 1:K2] <- t(ztw) %*% Z
## 	# XZ and ZX blocks
##   if (!block.diag) {
## 	  ztw2 <- 0*Z
## 	  for (k in 1:K2) ztw2[,k] <- Z[,k]*ret$h[,3]
## 	  h[K1 + 1:K2, 1:K1] <- t(ztw2)%*%X
## 	  h[1:K1, K1 + 1:K2] <- t(h[K1 + 1:K2, 1:K1])
## 	}
## 	
##   return (list(f = f, g = g, h = h))
## }


###################################################
### code chunk number 4: RegressionFactory.Rnw:350-351
###################################################
library(RegressionFactory)


###################################################
### code chunk number 5: RegressionFactory.Rnw:354-357
###################################################
loglike.logistic <- function(beta, X, y, fgh) {
  regfac.expand.1par(beta, X, y, fbase1.binomial.logit, fgh, n=1)
}


###################################################
### code chunk number 6: RegressionFactory.Rnw:360-368
###################################################
logprior.logistic <- function(beta, mu.beta, sd.beta, fgh) {
  f <- sum(dnorm(beta, mu.beta, sd.beta, log=TRUE))
  if (fgh==0) return (f)
  g <- -(beta-mu.beta)/sd.beta^2
  if (fgh==1) return (list(f=f, g=g))
  h <- diag(-1/sd.beta^2, nrow=length(beta))
  return (list(f=f, g=g, h=h))
}


###################################################
### code chunk number 7: RegressionFactory.Rnw:371-376
###################################################
logpost.logistic <- function(beta, X, y, mu.beta, sd.beta, fgh) {
  ret.loglike <- loglike.logistic(beta, X, y, fgh)
  ret.logprior <- logprior.logistic(beta, mu.beta, sd.beta, fgh)
  regfac.merge(ret.loglike, ret.logprior, fgh=fgh)
}


###################################################
### code chunk number 8: RegressionFactory.Rnw:381-387
###################################################
N <- 1000
K <- 5
X <- matrix(runif(N*K, min=-0.5, max=+0.5), ncol=K)
beta <- runif(K, min=-0.5, max=+0.5)
y <- rbinom(N, size = 1, prob = 1/(1+exp(-X%*%beta)))
beta.glm <- glm(y~X-1, family="binomial")$coefficients


###################################################
### code chunk number 9: RegressionFactory.Rnw:390-403
###################################################
library(sns)
nsmp <- 10
mu.beta <- 0.0
sd.beta <- 1000
beta.smp <- array(NA, dim=c(nsmp,K)) 
beta.tmp <- rep(0,K)
for (n in 1:nsmp) {
  beta.tmp <- sns(beta.tmp, fghEval=logpost.logistic, X=X, y=y
    , mu.beta=mu.beta, sd.beta=sd.beta, fgh=2, rnd=FALSE)
  beta.smp[n,] <- beta.tmp
}
beta.sns <- colMeans(beta.smp[(nsmp/2+1):nsmp,])
cbind(beta.glm, beta.sns)


###################################################
### code chunk number 10: RegressionFactory.Rnw:406-419
###################################################
J <- 20
mu.beta.hb <- runif(K, min=-0.5, max=+0.5)
sd.beta.hb <- runif(K, min=0.5, max=1.0)
X.hb <- list()
y.hb <- list()
beta.hb <- array(NA, dim=c(J,K))
for (k in 1:K) {
  beta.hb[,k] <- rnorm(J, mu.beta.hb[k], sd.beta.hb[k])
}
for (j in 1:J) {
  X.hb[[j]] <- matrix(runif(N*K, min=-0.5, max=+0.5), ncol=K)
  y.hb[[j]] <- rbinom(N, size=1, prob=1/(1+exp(-X%*%beta.hb[j,])))
}


###################################################
### code chunk number 11: RegressionFactory.Rnw:422-427
###################################################
beta.glm.all <- array(NA, dim=c(J,K))
for (j in 1:J) {
  beta.glm.all[j,] <- glm(y.hb[[j]]~X.hb[[j]]-1
    , family="binomial")$coefficients
}


###################################################
### code chunk number 12: RegressionFactory.Rnw:430-441
###################################################
beta.smp.hb <- array(NA, dim=c(nsmp,J,K)) 
beta.tmp.hb <- array(0.0, dim=c(J,K))
for (n in 1:nsmp) {
  for (j in 1:J) {
    beta.tmp.hb[j,] <- sns(beta.tmp.hb[j,], fghEval=logpost.logistic
      , X=X.hb[[j]], y=y.hb[[j]]
      , mu.beta=mu.beta.hb, sd.beta=sd.beta.hb, fgh=2, rnd=F)
  }
  beta.smp.hb[n,,] <- beta.tmp.hb
}
beta.sns.hb <- apply(beta.smp.hb[(nsmp/2+1):nsmp,,], c(2,3), mean)


###################################################
### code chunk number 13: RegressionFactory.Rnw:444-446
###################################################
head(beta.glm.all)
head(beta.sns.hb)


###################################################
### code chunk number 14: shrinkage_plot
###################################################
plot(beta.glm.all[,1], beta.sns.hb[,1]
     , xlab="Unpooled Coefficients"
     , ylab="Pooled Coefficients")
abline(a=0, b=1)


###################################################
### code chunk number 15: fig1
###################################################
plot(beta.glm.all[,1], beta.sns.hb[,1]
     , xlab="Unpooled Coefficients"
     , ylab="Pooled Coefficients")
abline(a=0, b=1)


###################################################
### code chunk number 16: RegressionFactory.Rnw:470-471
###################################################
library(dglm)


###################################################
### code chunk number 17: RegressionFactory.Rnw:474-480
###################################################
loglike.linreg <- function(coeff, X, y, fgh, vd = F) {
  if (vd) regfac.expand.2par(coeff = coeff, X = X, Z = X, y = y
    , fbase2 = fbase2.gaussian.identity.log, fgh = fgh, block.diag = F)
  else regfac.expand.2par(coeff = coeff, X = X, y = y
    , fbase2 = fbase2.gaussian.identity.log, fgh = fgh, block.diag = F)
}


###################################################
### code chunk number 18: RegressionFactory.Rnw:483-491
###################################################
N <- 1000
K <- 5
X <- matrix(runif(N*K, min=-0.5, max=+0.5), ncol=K)
beta <- runif(K, min=-0.5, max=+0.5)
gamma <- runif(K, min=-0.5, max=+0.5)
mean.vec <- X%*%beta
sd.vec <- exp(X%*%gamma)
y <- rnorm(N, mean.vec, sd.vec)


###################################################
### code chunk number 19: RegressionFactory.Rnw:494-502
###################################################
# constant-dispersion model
est.glm <- lm(y~X-1)
beta.glm <- est.glm$coefficients
sigma.glm <- summary(est.glm)$sigma
# varying-dispersion model
est.dglm <- dglm(y~X-1, dformula = ~X-1, family = "gaussian", dlink = "log")
beta.dglm <- est.dglm$coefficients
gamma.dglm <- est.dglm$dispersion.fit$coefficients


###################################################
### code chunk number 20: RegressionFactory.Rnw:505-529
###################################################
# constant-dispersion
coeff.smp <- array(NA, dim=c(nsmp, K+1)) 
coeff.tmp <- rep(0, K+1)
for (n in 1:nsmp) {
  coeff.tmp <- sns(coeff.tmp, fghEval=loglike.linreg
    , X=X, y=y, fgh=2, vd = F, rnd = F)
  coeff.smp[n,] <- coeff.tmp
}
beta.sns.cd <- colMeans(coeff.smp[(nsmp/2+1):nsmp, 1:K])
sigma.sns.cd <- sqrt(exp(mean(coeff.smp[(nsmp/2+1):nsmp, K+1])))
cbind(beta.glm, beta.sns.cd)
cbind(sigma.glm, sigma.sns.cd)
# varying-dispersion
coeff.smp <- array(NA, dim=c(nsmp, 2*K)) 
coeff.tmp <- rep(0, 2*K)
for (n in 1:nsmp) {
  coeff.tmp <- sns(coeff.tmp, fghEval=loglike.linreg
    , X=X, y=y, fgh=2, vd = T, rnd = F)
  coeff.smp[n,] <- coeff.tmp
}
beta.sns.vd <- colMeans(coeff.smp[(nsmp/2+1):nsmp, 1:K])
gamma.sns.vd <- colMeans(coeff.smp[(nsmp/2+1):nsmp, K+1:K])
cbind(beta.dglm, beta.sns.vd)
cbind(gamma.dglm, gamma.sns.vd)


###################################################
### code chunk number 21: RegressionFactory.Rnw:547-552
###################################################
N <- 1000
K <- 5
X <- matrix(runif(N*K, min=-0.5, max=+0.5), ncol=K)
beta <- runif(K, min=-0.5, max=+0.5)
y <- rgeom(N, prob = 1/(1+exp(-X%*%beta)))


###################################################
### code chunk number 22: RegressionFactory.Rnw:555-564
###################################################
loglike.geometric <- function(beta, X, y, fgh) {
  regfac.expand.1par(beta, X, y, fbase1.geometric.logit, fgh)
}
beta.est <- rep(0,K)
for (n in 1:10) {
  beta.est <- sns(beta.est, fghEval=loglike.geometric
    , X=X, y=y, fgh=2, rnd = F)
}
cbind(beta, beta.est)


