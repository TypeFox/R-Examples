### R code from vignette source 'MfUSampler.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: MfUSampler.Rnw:94-95
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: MfUSampler.Rnw:161-164
###################################################
library("MfUSampler")
my.seed <- 0
nsmp <- 10


###################################################
### code chunk number 3: MfUSampler.Rnw:172-177
###################################################
z <- c(1, 4, 7, 10, 13, 16, 19, 24)
m1.prior <- c(17, 26, 39, 27, 35, 37, 26, 23)
m2.prior <- c(215, 218, 137, 62, 36, 16, 13, 15)
m1.current <- c(46, 52, 44, 54, 38, 39, 23, 52)
m2.current <- c(290, 211, 134, 91, 53, 42, 23, 32)


###################################################
### code chunk number 4: MfUSampler.Rnw:180-181
###################################################
X <- cbind(1, z, z^2)


###################################################
### code chunk number 5: MfUSampler.Rnw:187-199
###################################################
loglike <- function(beta, X, m1, m2) {
  beta <- as.numeric(beta)
  Xbeta <- X %*% beta
  return (-sum((m1 + m2) * log(1 + exp(-Xbeta)) + m2 * Xbeta))
}
logprior <- function(beta, beta0 , W) {
  return (-0.5 * t(beta - beta0) %*% solve(W) %*% (beta - beta0))
}
logpost <- function(beta, X, m1, m2
  , beta0 = rep(0,0, 3), W = diag(1e+6, nrow = 3)) {
  return (logprior(beta, beta0, W) + loglike(beta, X, m1, m2))
}


###################################################
### code chunk number 6: MfUSampler.Rnw:202-204
###################################################
m1.total <- m1.prior + m1.current
m2.total <- m2.prior + m2.current


###################################################
### code chunk number 7: MfUSampler.Rnw:207-212 (eval = FALSE)
###################################################
## set.seed(my.seed)
## beta.ini <- c(0.0, 0.0, 0.0)
## beta.smp <- MfU.Sample.Run(beta.ini, logpost, nsmp = nsmp
##   , X = X, m1 = m1.total, m2 = m2.total)
## summ.slice <- summary(beta.smp)


###################################################
### code chunk number 8: MfUSampler.Rnw:214-215
###################################################
load("summ.slice")


###################################################
### code chunk number 9: MfUSampler.Rnw:217-218
###################################################
print(summ.slice)


###################################################
### code chunk number 10: MfUSampler.Rnw:221-222
###################################################
print(summ.slice$covar)


###################################################
### code chunk number 11: MfUSampler.Rnw:227-244
###################################################
logpost.fg <- function(beta, X, m1, m2
  , beta0 = rep(0.0, 3), W = diag(1e+3, nrow = 3)
  , grad = FALSE) {
  Xbeta <- X %*% beta
  
  if (grad) {
    log.prior.d <- -solve(W) %*% (beta - beta0)
    log.like.d <- t(X) %*% ((m1 + m2) / (1 + exp(Xbeta)) - m2)
    return (log.prior.d + log.like.d)
  }
  
  log.prior <- -0.5 * t(beta - beta0) %*% solve(W) %*% (beta - beta0)
  log.like <- -sum((m1 + m2) * log(1 + exp(-Xbeta)) + m2 * Xbeta)
  log.post <- log.prior + log.like

  return (log.post)
}


###################################################
### code chunk number 12: MfUSampler.Rnw:247-255
###################################################
set.seed(my.seed)
beta.ini <- c(0.0, 0.0, 0.0)
beta.smp <- MfU.Sample.Run(beta.ini, logpost.fg, nsmp = nsmp
  , uni.sampler = "ars"
  , control = MfU.Control(3, ars.x = list(c(-10, 0, 10)
                , c(-1, 0, 1), c(-0.1, 0.0, 0.1)))
  , X = X, m1 = m1.total, m2 = m2.total)
summ.ars <- summary(beta.smp)


###################################################
### code chunk number 13: MfUSampler.Rnw:257-258
###################################################
load("summ.ars")


###################################################
### code chunk number 14: MfUSampler.Rnw:260-261
###################################################
print(summ.ars)


###################################################
### code chunk number 15: MfUSampler.Rnw:267-270
###################################################
predfunc.mean <- function(beta, X) {
  return (1/(1 + exp(-X %*% beta)))
}


###################################################
### code chunk number 16: MfUSampler.Rnw:273-276
###################################################
pred.mean <- predict(beta.smp, predfunc.mean, X)
predmean.summ <- summary(pred.mean)
print(predmean.summ, n = 8)


###################################################
### code chunk number 17: MfUSampler.Rnw:279-285
###################################################
predfunc.binary <- function(beta, X) {
  return (1*(runif(nrow(X)) < 1/(1 + exp(-X %*% beta))))
}
pred.binary <- predict(beta.smp, predfunc.binary, X)
predbinary.summ <- summary(pred.binary)
print(predbinary.summ, n = 8)


###################################################
### code chunk number 18: MfUSampler.Rnw:294-310
###################################################
library("mvtnorm")
set.seed(my.seed)
nrep <- 50
m.current <- m1.current + m2.current
nz.exp <- nrep * sum(m.current)
jitter <- 1.0
z.exp <- sample(z, size = nz.exp, replace = T, prob = m.current) +
  (2*runif(nz.exp) - 1) * jitter
X.exp <- cbind(1, z.exp, z.exp^2)
beta0.prior <- c(-3.17, 0.33, -0.007)
W.prior <- 1e-4 * matrix(c(638, -111, 3.9
  , -111, 24.1, -0.9, 3.9, -0.9, 0.04), ncol = 3)
ngrp <- 50
beta.mat <- t(rmvnorm(ngrp, mean = beta0.prior, sigma = W.prior))
y.mat.exp <- 1* matrix(runif(ngrp * nz.exp) < 
  1 / (1 + exp(-X.exp %*% beta.mat)), ncol = ngrp)


###################################################
### code chunk number 19: MfUSampler.Rnw:316-332
###################################################
hb.logprior <- function(beta.flat, beta0, W) {
  beta.mat <- matrix(beta.flat, nrow = 3)
  return (sum(apply(beta.mat, 2, logprior, beta0, W)))
}
hb.loglike <- function(beta.flat, X, y) {
  beta.mat <- matrix(beta.flat, nrow = 3)
  ngrp <- ncol(beta.mat)
  return (sum(sapply(1:ngrp, function(n) {
    xbeta <- X %*% beta.mat[, n]
    return (-sum((1-y[, n]) * xbeta + log(1 + exp(-xbeta))))
  })))
}
hb.logpost <- function(beta.flat, X, y, beta0, W) {
  return (hb.logprior(beta.flat, beta0, W) + 
            hb.loglike(beta.flat, X, y))
}


###################################################
### code chunk number 20: MfUSampler.Rnw:335-342 (eval = FALSE)
###################################################
## nsmp <- 10
## set.seed(my.seed)
## beta.flat.ini <- rep(0.0, 3 * ngrp)
## beta.flat.smp <- MfU.Sample.Run(beta.flat.ini, hb.logpost
##                      , X = X.exp, y = y.mat.exp
##                      , beta0 = beta0.prior, W = W.prior, nsmp = nsmp)
## t.naive <- attr(beta.flat.smp, "t")


###################################################
### code chunk number 21: MfUSampler.Rnw:344-345
###################################################
load("t.naive")


###################################################
### code chunk number 22: MfUSampler.Rnw:347-348
###################################################
cat("hb sampling time - naive method:", t.naive, "sec\n")


###################################################
### code chunk number 23: MfUSampler.Rnw:353-364
###################################################
hb.loglike.grp <- function(beta, X, y) {
  beta <- as.numeric(beta)
  xbeta <- X %*% beta
  return (-sum((1-y) * xbeta + log(1 + exp(-xbeta))))
}
hb.logprior.grp <- logprior
hb.logpost.grp <- function(beta, X, y
  , beta0 = rep(0,0, 3), W = diag(1e+6, nrow = 3)) {
  return (hb.logprior.grp(beta, beta0, W) +
            hb.loglike.grp(beta, X, y))
}


###################################################
### code chunk number 24: MfUSampler.Rnw:367-368
###################################################
t.revised <- 10.7


###################################################
### code chunk number 25: MfUSampler.Rnw:370-383 (eval = FALSE)
###################################################
## set.seed(my.seed)
## beta.mat.buff <- matrix(rep(0.0, 3 * ngrp), nrow = 3)
## beta.mat.smp <- array(NA, dim = c(nsmp, 3, ngrp))
## t.revised <- proc.time()[3]
## for (i in 1:nsmp) {
##   for (n in 1:ngrp) {
##     beta.mat.buff[, n] <- MfU.Sample(beta.mat.buff[, n], hb.logpost.grp
##       , uni.sampler = "slice", X = X.exp
##       , y = y.mat.exp[, n], beta0 = beta0.prior, W = W.prior)
##   }
##   beta.mat.smp[i, , ] <- beta.mat.buff
## }
## t.revised <- proc.time()[3] - t.revised


###################################################
### code chunk number 26: MfUSampler.Rnw:385-386
###################################################
load("t.revised")


###################################################
### code chunk number 27: MfUSampler.Rnw:388-389
###################################################
cat("hb sampling time - revised method:", t.revised, "sec\n")


###################################################
### code chunk number 28: MfUSampler.Rnw:392-410 (eval = FALSE)
###################################################
## library("doParallel")
## ncores <- 2
## registerDoParallel(ncores)
## 
## set.seed(my.seed)
## beta.mat.buff <- matrix(rep(0.0, 3 * ngrp), nrow = 3)
## beta.mat.smp <- array(NA, dim = c(nsmp, 3, ngrp))
## t.parallel <- proc.time()[3]
## for (i in 1:nsmp) {
##   beta.mat.buff <- foreach(n=1:ngrp, .combine = cbind
##     , .options.multicore=list(preschedule=TRUE)) %dopar% {
##       MfU.Sample(beta.mat.buff[, n], hb.logpost.grp, uni.sampler = "slice"
##         , X = X.exp, y = y.mat.exp[, n]
##         , beta0 = beta0.prior, W = W.prior)
##   }
##   beta.mat.smp[i, , ] <- beta.mat.buff
## }
## t.parallel <- proc.time()[3] - t.parallel


###################################################
### code chunk number 29: MfUSampler.Rnw:412-413
###################################################
load("t.parallel")


###################################################
### code chunk number 30: MfUSampler.Rnw:415-416
###################################################
cat("hb sampling time - revised & parallel method:", t.parallel, "sec\n")


###################################################
### code chunk number 31: MfUSampler.Rnw:422-445
###################################################
library("RcppArmadillo")
library("inline")
code <- "
  arma::vec beta_cpp = Rcpp::as<arma::vec>(beta);
  arma::mat X_cpp = Rcpp::as<arma::mat>(X);
  arma::vec y_cpp = Rcpp::as<arma::vec>(y);
  arma::vec xbeta = X_cpp * beta_cpp;
  int n = X_cpp.n_rows;
  double logp = 0.0;
  for (int i=0; i<n; i++) {
    // (1-y[, n]) * xbeta + log(1 + exp(-xbeta))
    logp -= (1.0 - y_cpp[i]) * xbeta[i] + log(1.0 + exp(-xbeta[i]));
  }
  return Rcpp::wrap(logp);
"
hb.loglike.grp.rcpp <- cxxfunction(
  signature(beta = "numeric", X = "numeric", y = "numeric")
  , code, plugin="RcppArmadillo")
hb.logpost.grp.rcpp <- function(beta, X, y
  , beta0 = rep(0,0, 3), W = diag(1e+6, nrow = 3)) {
  return (hb.logprior.grp(beta, beta0, W) +
            hb.loglike.grp.rcpp(beta, X, y))
}


###################################################
### code chunk number 32: MfUSampler.Rnw:448-462 (eval = FALSE)
###################################################
## set.seed(my.seed)
## beta.mat.buff <- matrix(rep(0.0, 3 * ngrp), nrow = 3)
## beta.mat.smp <- array(NA, dim = c(nsmp, 3, ngrp))
## t.rcpp <- proc.time()[3]
## for (i in 1:nsmp) {
##   beta.mat.buff <- foreach(n=1:ngrp
##     , .combine = cbind, .options.multicore=list(preschedule=TRUE)) %dopar% {
##       MfU.Sample(beta.mat.buff[, n], hb.logpost.grp.rcpp, uni.sampler = "slice"
##         , X = X.exp, y = y.mat.exp[, n]
##         , beta0 = beta0.prior, W = W.prior)
##   }
##   beta.mat.smp[i, , ] <- beta.mat.buff
## }
## t.rcpp <- proc.time()[3] - t.rcpp


###################################################
### code chunk number 33: MfUSampler.Rnw:464-465
###################################################
load("t.rcpp")


###################################################
### code chunk number 34: MfUSampler.Rnw:467-468
###################################################
cat("hb sampling time - revised & parallel & rcpp method:", t.rcpp, "sec\n")


###################################################
### code chunk number 35: my_plot
###################################################
yplot <- c(t.naive, t.revised, t.parallel, t.rcpp)
names(yplot) <- c("naive", " + separability", " + parallel", " + Rcpp")
xmids <- barplot(yplot, xlab = "optimization level", yaxt = "n", space = 0.5
                 #, log = "y"
                 #, ylab = "time (sec)"
                 , ylim = c(0, yplot[1] + 10)
        )
text(x = xmids, y = yplot/2, labels = format(yplot, digits = 3))
text(x = xmids[2:4], y = yplot[2:4] + 8, labels = paste(format(c(t.naive / t.revised, t.revised / t.parallel, t.parallel / t.rcpp), digits = 2), "x", sep = ""))


###################################################
### code chunk number 36: fig1
###################################################
yplot <- c(t.naive, t.revised, t.parallel, t.rcpp)
names(yplot) <- c("naive", " + separability", " + parallel", " + Rcpp")
xmids <- barplot(yplot, xlab = "optimization level", yaxt = "n", space = 0.5
                 #, log = "y"
                 #, ylab = "time (sec)"
                 , ylim = c(0, yplot[1] + 10)
        )
text(x = xmids, y = yplot/2, labels = format(yplot, digits = 3))
text(x = xmids[2:4], y = yplot[2:4] + 8, labels = paste(format(c(t.naive / t.revised, t.revised / t.parallel, t.parallel / t.rcpp), digits = 2), "x", sep = ""))


###################################################
### code chunk number 37: MfUSampler.Rnw:504-505
###################################################
sessionInfo()


