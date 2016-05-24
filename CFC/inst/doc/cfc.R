### R code from vignette source 'cfc.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: cfc.Rnw:96-97
###################################################
options(prompt = "R> ", continue = "+  ", width = 80, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: firstplot
###################################################
seed.no <- 0
set.seed(seed.no)

old.par <- par(mfrow = c(1,2))

library("BSGW")
library("Hmisc")
library("CFC")
data(bmt)

predict.S.core <- function(t, xbeta, xbetas, alpha.min, alpha.max, ordweib) {
  if (t == 0.0) return (1.0)
  if (ordweib) {
    alpha <- exp(xbetas)
  } else {
    alpha <- alpha.min + (alpha.max-alpha.min)/(1+exp(-xbetas))
  }
  H <- t^alpha * exp(xbeta)
  S <- exp(-H)
  return (S)
}

dataset <- "ovarian"

operator <- mean

iter <- 5000
ordweib <- TRUE

lung <- lung[complete.cases(lung), ]

est <- bsgw(Surv(futime, fustat) ~ as.factor(ecog.ps) + as.factor(rx), ovarian, ordweib = ordweib
            , control=bsgw.control(iter=iter, nskip=100), print.level = 0)
est.ref <- survreg(Surv(futime, fustat) ~ as.factor(ecog.ps) + as.factor(rx), ovarian)

X <- est$X
Xs <- est$Xs
beta.smp <- est$smp$beta
betas.smp <- est$smp$betas
alpha.min <- est$control$alpha.min
alpha.max <- est$control$alpha.max

beta.avg <- apply(beta.smp[(iter/2+1):iter, ], 2, operator)
betas.avg <- apply(betas.smp[(iter/2+1):iter, , drop = F], 2, operator)

t.compare <- 0.5 * est$tmax

xbeta.smp <- X %*% t(beta.smp)
xbetas.smp <- Xs %*% t(betas.smp)

xbetas.avg <- log(operator(exp(xbetas.smp[2, ])))

S.wrong <- predict.S.core(t.compare, X %*% beta.avg, xbetas.avg, alpha.min, alpha.max, ordweib)

S.smp <- sapply(1:iter, function(i) {
  predict.S.core(t.compare, xbeta.smp[, i], xbetas.smp[, i], alpha.min, alpha.max, ordweib)
})
S.right <- apply(S.smp[, (iter/2+1):iter], 1, operator)
S.right.sd <- apply(S.smp[, (iter/2+1):iter], 1, sd)
S.right.se <- S.right.sd / sqrt(iter/2)
S.right.plus <- S.right + 2 * S.right.se
S.right.minus <- S.right - 2 * S.right.se

my.range.x <- range(S.wrong, S.right)
my.range.y <- my.range.x

errbar(x = S.wrong, y = S.right, yplus = S.right.plus, yminus = S.right.minus
       , ylim = my.range.y, xlim = my.range.x
       , ylab = "Bayesian survival probability", xlab = "non-Bayesian survival probability"
       )
abline(a = 0, b = 1)

library("survival")
days.in.year <- 365.25
fitKM <- survfit(Surv(stop, event=='pcm') ~1, data=mgus1,
                    subset=(start==0))
fitCI <- survfit(Surv(stop, status*as.numeric(event), type="mstate") ~1,
                    data=mgus1, subset=(start==0))

plot(fitCI$time / days.in.year, fitCI$prev[,1], type = "l"
     , ylim = c(0.0, 0.4), xlim = c(0.0, 7300) / days.in.year
     , xlab = "Years post diagnosis of MGUS", ylab = "Incidence Probability"
     , lwd = 2)
lines(fitKM$time / days.in.year, 1 - fitKM$surv, lty = 2, lwd = 2)
legend("topleft", legend = c("Competing-Risk", "Kaplan-Meyer"), lty = 1:2, lwd = rep(2, 1))

par(old.par)


###################################################
### code chunk number 3: twofigs
###################################################
seed.no <- 0
set.seed(seed.no)

old.par <- par(mfrow = c(1,2))

library("BSGW")
library("Hmisc")
library("CFC")
data(bmt)

predict.S.core <- function(t, xbeta, xbetas, alpha.min, alpha.max, ordweib) {
  if (t == 0.0) return (1.0)
  if (ordweib) {
    alpha <- exp(xbetas)
  } else {
    alpha <- alpha.min + (alpha.max-alpha.min)/(1+exp(-xbetas))
  }
  H <- t^alpha * exp(xbeta)
  S <- exp(-H)
  return (S)
}

dataset <- "ovarian"

operator <- mean

iter <- 5000
ordweib <- TRUE

lung <- lung[complete.cases(lung), ]

est <- bsgw(Surv(futime, fustat) ~ as.factor(ecog.ps) + as.factor(rx), ovarian, ordweib = ordweib
            , control=bsgw.control(iter=iter, nskip=100), print.level = 0)
est.ref <- survreg(Surv(futime, fustat) ~ as.factor(ecog.ps) + as.factor(rx), ovarian)

X <- est$X
Xs <- est$Xs
beta.smp <- est$smp$beta
betas.smp <- est$smp$betas
alpha.min <- est$control$alpha.min
alpha.max <- est$control$alpha.max

beta.avg <- apply(beta.smp[(iter/2+1):iter, ], 2, operator)
betas.avg <- apply(betas.smp[(iter/2+1):iter, , drop = F], 2, operator)

t.compare <- 0.5 * est$tmax

xbeta.smp <- X %*% t(beta.smp)
xbetas.smp <- Xs %*% t(betas.smp)

xbetas.avg <- log(operator(exp(xbetas.smp[2, ])))

S.wrong <- predict.S.core(t.compare, X %*% beta.avg, xbetas.avg, alpha.min, alpha.max, ordweib)

S.smp <- sapply(1:iter, function(i) {
  predict.S.core(t.compare, xbeta.smp[, i], xbetas.smp[, i], alpha.min, alpha.max, ordweib)
})
S.right <- apply(S.smp[, (iter/2+1):iter], 1, operator)
S.right.sd <- apply(S.smp[, (iter/2+1):iter], 1, sd)
S.right.se <- S.right.sd / sqrt(iter/2)
S.right.plus <- S.right + 2 * S.right.se
S.right.minus <- S.right - 2 * S.right.se

my.range.x <- range(S.wrong, S.right)
my.range.y <- my.range.x

errbar(x = S.wrong, y = S.right, yplus = S.right.plus, yminus = S.right.minus
       , ylim = my.range.y, xlim = my.range.x
       , ylab = "Bayesian survival probability", xlab = "non-Bayesian survival probability"
       )
abline(a = 0, b = 1)

library("survival")
days.in.year <- 365.25
fitKM <- survfit(Surv(stop, event=='pcm') ~1, data=mgus1,
                    subset=(start==0))
fitCI <- survfit(Surv(stop, status*as.numeric(event), type="mstate") ~1,
                    data=mgus1, subset=(start==0))

plot(fitCI$time / days.in.year, fitCI$prev[,1], type = "l"
     , ylim = c(0.0, 0.4), xlim = c(0.0, 7300) / days.in.year
     , xlab = "Years post diagnosis of MGUS", ylab = "Incidence Probability"
     , lwd = 2)
lines(fitKM$time / days.in.year, 1 - fitKM$surv, lty = 2, lwd = 2)
legend("topleft", legend = c("Competing-Risk", "Kaplan-Meyer"), lty = 1:2, lwd = rep(2, 1))

par(old.par)


###################################################
### code chunk number 4: cfc.Rnw:376-383
###################################################
library("CFC")
data(bmt)
rel.tol <- 1e-3
idx.train <- sample(1:nrow(bmt), size = 0.7 * nrow(bmt))
idx.pred <- setdiff(1:nrow(bmt), idx.train)
nobs.train <- length(idx.train)
nobs.pred <- length(idx.pred)


###################################################
### code chunk number 5: cfc.Rnw:386-388
###################################################
out.weib <- cfc.survreg(Surv(time, cause) ~ platelet + age + tcell, 
  bmt[idx.train, ], bmt[idx.pred, ], rel.tol = rel.tol)


###################################################
### code chunk number 6: cfc.Rnw:391-392
###################################################
summ <- summary(out.weib, obs.idx = which(bmt$age[idx.pred] > 0))


###################################################
### code chunk number 7: myplot
###################################################
plot(summ, which = 1)


###################################################
### code chunk number 8: fig1
###################################################
plot(summ, which = 1)


###################################################
### code chunk number 9: myplot2
###################################################
old.par <- par(mfrow=c(1,2)); plot(summ, which = 2); par(old.par)


###################################################
### code chunk number 10: fig2
###################################################
old.par <- par(mfrow=c(1,2)); plot(summ, which = 2); par(old.par)


###################################################
### code chunk number 11: cfc.Rnw:420-423
###################################################
out.expo <- cfc.survreg(Surv(time, cause) ~ platelet + age + tcell, 
  bmt[idx.train, ], bmt[idx.pred, ],
  dist = "exponential", rel.tol = rel.tol)


###################################################
### code chunk number 12: cfc.Rnw:426-429
###################################################
out.mix <- cfc.survreg(Surv(time, cause) ~ platelet + age + tcell, 
  bmt[idx.train, ], bmt[idx.pred, ],
  dist = c("weibull", "exponential"), rel.tol = rel.tol)


###################################################
### code chunk number 13: cfc.Rnw:436-441
###################################################
out.prep <- cfc.prepdata(Surv(time, cause) ~ platelet + age + tcell, bmt)
f1 <- out.prep$formula.list[[1]]
f2 <- out.prep$formula.list[[2]]
dat <- out.prep$dat
tmax <- out.prep$tmax


###################################################
### code chunk number 14: cfc.Rnw:444-452
###################################################
library("BSGW")
nsmp <- 10
reg1 <- bsgw(f1, dat[idx.train, ],
  control = bsgw.control(iter = nsmp),
  ordweib = T, print.level = 0)
reg2 <- bsgw(f2, dat[idx.train, ], 
  control = bsgw.control(iter = nsmp),
  ordweib = T, print.level = 0)


###################################################
### code chunk number 15: cfc.Rnw:455-473
###################################################
X.pred <- as.matrix(cbind(1, bmt[idx.pred, c("platelet", "age", "tcell")]))
arg.1 <- list(nobs = nobs.pred, natt = 4, nsmp = nsmp,
  X = X.pred, alpha = exp(reg1$smp$betas),
  beta = reg1$smp$beta)
arg.2 <- list(nobs = nobs.pred, natt = 4, nsmp = nsmp, 
  X = X.pred, alpha = exp(reg2$smp$betas), 
  beta = reg2$smp$beta)
arg.list <- list(arg.1, arg.2)
survfunc <- function(t, args, n) {
  nobs <- args$nobs; natt <- args$natt; nsmp <- args$nsmp
  alpha <- args$alpha; beta <- args$beta; X <- args$X
  idx.smp <- floor((n - 1) / nobs) + 1
  idx.obs <- n - (idx.smp - 1) * nobs
  return (exp(- t ^ alpha[idx.smp] * 
                exp(sum(X[idx.obs, ] * 
                          beta[idx.smp, ]))));
}
f.list <- list(survfunc, survfunc)


###################################################
### code chunk number 16: cfc.Rnw:476-483
###################################################
rel.tol <- 1e-4
tout <- seq(from = 0.0, to = tmax, length.out = 10)
t.R <- proc.time()[3]
out.cfc.R <- cfc(f.list, arg.list, nobs.pred * nsmp, tout,
  rel.tol = rel.tol)
t.R <- proc.time()[3] - t.R
cat("t.R:", t.R, "sec\n")


###################################################
### code chunk number 17: cfc.Rnw:487-494
###################################################
ncores <- 2
tout <- seq(from = 0.0, to = tmax, length.out = 10)
t.R.par <- proc.time()[3]
out.cfc.R.par <- cfc(f.list, arg.list, nobs.pred * nsmp, tout,
  rel.tol = rel.tol, ncores = ncores)
t.R.par <- proc.time()[3] - t.R.par
cat("t.R.par:", t.R.par, "sec\n")


###################################################
### code chunk number 18: cfc.Rnw:497-498
###################################################
cat("parallelization speedup - R:", t.R / t.R.par, "\n")


###################################################
### code chunk number 19: cfc.Rnw:566-577
###################################################
tout <- seq(from = 0.0, to = tmax, length.out = 10)
library("Rcpp")
Rcpp::sourceCpp("weib.cpp")
f.list.Cpp.1 <- list(weib_getPtr_func(), weib_getPtr_init(),
  weib_getPtr_free())
f.list.Cpp <- list(f.list.Cpp.1, f.list.Cpp.1)
t.Cpp <- proc.time()[3]
out.cfc.Cpp <- cfc(f.list.Cpp, arg.list, nobs.pred * nsmp, tout,
  rel.tol = rel.tol)
t.Cpp <- proc.time()[3] - t.Cpp
cat("t.Cpp:", t.Cpp, "sec\n")


###################################################
### code chunk number 20: cfc.Rnw:580-581
###################################################
all.equal(out.cfc.R, out.cfc.Cpp)


###################################################
### code chunk number 21: cfc.Rnw:584-585
###################################################
cat("C++-vs-R speedup:", t.R / t.Cpp, "\n")


###################################################
### code chunk number 22: cfc.Rnw:588-607
###################################################
nsmp <- 1000
reg1 <- bsgw(f1, dat[idx.train, ],
  control = bsgw.control(iter = nsmp),
  ordweib = T, print.level = 0)
reg2 <- bsgw(f2, dat[idx.train, ], 
  control = bsgw.control(iter = nsmp),
  ordweib = T, print.level = 0)
arg.1 <- list(nobs = nobs.pred, natt = 4, nsmp = nsmp,
  X = X.pred, alpha = exp(reg1$smp$betas),
  beta = reg1$smp$beta)
arg.2 <- list(nobs = nobs.pred, natt = 4, nsmp = nsmp, 
  X = X.pred, alpha = exp(reg2$smp$betas), 
  beta = reg2$smp$beta)
arg.list <- list(arg.1, arg.2)
t.Cpp <- proc.time()[3]
out.cfc.Cpp <- cfc(f.list.Cpp, arg.list, nobs.pred * nsmp, tout,
  rel.tol = rel.tol)
t.Cpp <- proc.time()[3] - t.Cpp
cat("t.Cpp:", t.Cpp, "sec\n")


###################################################
### code chunk number 23: cfc.Rnw:610-615
###################################################
t.Cpp.par <- proc.time()[3]
out.cfc.Cpp <- cfc(f.list.Cpp, arg.list, nobs.pred * nsmp, tout,
  rel.tol = rel.tol, ncores = ncores)
t.Cpp.par <- proc.time()[3] - t.Cpp.par
cat("t.Cpp.par:", t.Cpp.par, "sec\n")


###################################################
### code chunk number 24: cfc.Rnw:618-619
###################################################
cat("parallelization speedup - C++:", t.Cpp / t.Cpp.par, "\n")


###################################################
### code chunk number 25: cfc.Rnw:629-634
###################################################
prep <- cfc.prepdata(Surv(time, cause) ~ platelet + age + tcell, bmt)
f1 <- prep$formula.list[[1]]
f2 <- prep$formula.list[[2]]
dat <- prep$dat
tmax <- prep$tmax


###################################################
### code chunk number 26: cfc.Rnw:637-639
###################################################
library("survival")
reg1 <- survreg(f1, dat, x = TRUE)


###################################################
### code chunk number 27: cfc.Rnw:642-650
###################################################
library("randomForestSRC")
reg2 <- rfsrc(f2, dat)
rfsrc.survfunc <- function(t, args, n) {
  which.zero <- which(t < .Machine$double.eps)
  ret <- approx(args$time.interest, args$survival[n, ], t, rule = 2)$y
  ret[which.zero] <- 1.0
  return (ret)
}


###################################################
### code chunk number 28: cfc.Rnw:653-657
###################################################
f.list <- list(cfc.survreg.survprob, rfsrc.survfunc)
arg.list <- list(reg1, reg2)
tout <- seq(0.0, tmax, length.out = 10)
cfc.out <- cfc(f.list, arg.list, nrow(bmt), tout, rel.tol = 1e-3)


###################################################
### code chunk number 29: cfc.Rnw:681-682
###################################################
sessionInfo()


