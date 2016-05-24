### R code from vignette source 'bfst.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: foo
###################################################
options(keep.source = TRUE, width = 65)


###################################################
### code chunk number 2: library
###################################################
library(mcmc)


###################################################
### code chunk number 3: baz
###################################################
baz <- library(help = "mcmc")
baz <- baz$info[[1]]
baz <- baz[grep("Version", baz)]
baz <- sub("^Version: *", "", baz)
bazzer <- paste(R.version$major, R.version$minor, sep = ".")


###################################################
### code chunk number 4: set-seed
###################################################
set.seed(42)


###################################################
### code chunk number 5: frequentist
###################################################
data(logit)
out <- glm(y ~ x1 + x2 + x3 + x4, data = logit,
    family = binomial, x = TRUE)
summary(out)


###################################################
### code chunk number 6: models
###################################################
varnam <- names(coefficients(out))
varnam <- varnam[varnam != "(Intercept)"]
nvar <- length(varnam)

models <- NULL
foo <- seq(0, 2^nvar - 1) 
for (i in 1:nvar) {
    bar <- foo %/% 2^(i - 1)
    bar <- bar %% 2
    models <- cbind(bar, models, deparse.level = 0)
}
colnames(models) <- varnam
models


###################################################
### code chunk number 7: neighbor
###################################################
neighbors <- matrix(FALSE, nrow(models), nrow(models))
for (i in 1:nrow(neighbors)) {
    for (j in 1:ncol(neighbors)) {
        foo <- models[i, ]
        bar <- models[j, ]
        if (sum(foo != bar) == 1) neighbors[i, j] <- TRUE
    }
}


###################################################
### code chunk number 8: ludfun
###################################################
modmat <- out$x
y <- logit$y

ludfun <- function(state, log.pseudo.prior) {
    stopifnot(is.numeric(state))
    stopifnot(length(state) == ncol(models) + 2)
    icomp <- state[1]
    stopifnot(icomp == as.integer(icomp))
    stopifnot(1 <= icomp && icomp <= nrow(models))
    stopifnot(is.numeric(log.pseudo.prior))
    stopifnot(length(log.pseudo.prior) == nrow(models))
    beta <- state[-1]
    inies <- c(TRUE, as.logical(models[icomp, ]))
    beta.logl <- beta
    beta.logl[! inies] <- 0
    eta <- as.numeric(modmat %*% beta.logl)
    logp <- ifelse(eta < 0, eta - log1p(exp(eta)), - log1p(exp(- eta)))
    logq <- ifelse(eta < 0, - log1p(exp(eta)), - eta - log1p(exp(- eta)))
    logl <- sum(logp[y == 1]) + sum(logq[y == 0])
    logl + sum(dnorm(beta, 0, 2, log = TRUE)) + log.pseudo.prior[icomp]
}


###################################################
### code chunk number 9: try1
###################################################
state.initial <- c(nrow(models), out$coefficients)

qux <- rep(0, nrow(models))

out <- temper(ludfun, initial = state.initial, neighbors = neighbors,
    nbatch = 1000, blen = 100, log.pseudo.prior = qux)

names(out)
out$time


###################################################
### code chunk number 10: what
###################################################
ibar <- colMeans(out$ibatch)
ibar


###################################################
### code chunk number 11: adjust
###################################################
qux <- qux + pmin(log(max(ibar) / ibar), 10)
qux <- qux - min(qux)
qux


###################################################
### code chunk number 12: iterate
###################################################
lout <- suppressWarnings(try(load("bfst1.rda"), silent = TRUE))
if (inherits(lout, "try-error")) {
    qux.save <- qux
    time.save <- out$time
    repeat{
        out <- temper(out, log.pseudo.prior = qux)
        ibar <- colMeans(out$ibatch)
        qux <- qux + pmin(log(max(ibar) / ibar), 10)
        qux <- qux - min(qux)
        qux.save <- rbind(qux.save, qux, deparse.level = 0)
        time.save <- rbind(time.save, out$time, deparse.level = 0)
        if (max(ibar) / min(ibar) < 2) break
    }
    save(out, qux, qux.save, time.save, file = "bfst1.rda")
} else {
    .Random.seed <- out$final.seed
}
print(qux.save, digits = 3)
print(qux, digits = 3)
apply(time.save, 2, sum)


###################################################
### code chunk number 13: accept-i-x
###################################################
print(out$accepti, digits = 3)
print(out$acceptx, digits = 3)


###################################################
### code chunk number 14: accept-i-min
###################################################
min(as.vector(out$accepti), na.rm = TRUE)


###################################################
### code chunk number 15: scale
###################################################
out <- temper(out, scale = 0.5, log.pseudo.prior = qux)
time.save <- rbind(time.save, out$time, deparse.level = 0)
print(out$acceptx, digits = 3)


###################################################
### code chunk number 16: try6
###################################################
lout <- suppressWarnings(try(load("bfst2.rda"), silent = TRUE))
if (inherits(lout, "try-error")) {
    out <- temper(out, blen = 10 * out$blen, log.pseudo.prior = qux)
    save(out, file = "bfst2.rda")
} else {
    .Random.seed <- out$final.seed
}
time.save <- rbind(time.save, out$time, deparse.level = 0)
foo <- apply(time.save, 2, sum)
foo.min <- floor(foo[1] / 60)
foo.sec <- foo[1] - 60 * foo.min
c(foo.min, foo.sec)


###################################################
### code chunk number 17: doit
###################################################
log.10.unnorm.bayes <- (qux - log(colMeans(out$ibatch))) / log(10)
k <- seq(along = log.10.unnorm.bayes)[log.10.unnorm.bayes
    == min(log.10.unnorm.bayes)]
models[k, ]

log.10.bayes <- log.10.unnorm.bayes - log.10.unnorm.bayes[k]
log.10.bayes


###################################################
### code chunk number 18: doit-se-one
###################################################
fred <- var(out$ibatch) / out$nbatch
sally <- colMeans(out$ibatch)
mcse.log.10.bayes <- (1 / log(10)) * sqrt(diag(fred) / sally^2 -
    2 * fred[ , k] / (sally * sally[k]) +
    fred[k, k] / sally[k]^2)
mcse.log.10.bayes

foompter <- cbind(models, log.10.bayes, mcse.log.10.bayes)
round(foompter, 5)


###################################################
### code chunk number 19: doit-too
###################################################
ibar <- colMeans(out$ibatch)
herman <- sweep(out$ibatch, 2, ibar, "/")
herman <- sweep(herman, 1, herman[ , k], "-")
mcse.log.10.bayes.too <- (1 / log(10)) *
    apply(herman, 2, sd) /sqrt(out$nbatch)
all.equal(mcse.log.10.bayes, mcse.log.10.bayes.too)


