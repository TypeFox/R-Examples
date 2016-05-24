### R code from vignette source 'demo.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: foo
###################################################
options(keep.source = TRUE, width = 60)
foo <- packageDescription("mcmc")


###################################################
### code chunk number 2: frequentist
###################################################
library(mcmc)
data(logit)
out <- glm(y ~ x1 + x2 + x3 + x4, data = logit,
    family = binomial(), x = TRUE)
summary(out)


###################################################
### code chunk number 3: log.unnormalized.posterior
###################################################
x <- out$x
y <- out$y

lupost <- function(beta, x, y) {
    eta <- as.numeric(x %*% beta)
    logp <- ifelse(eta < 0, eta - log1p(exp(eta)), - log1p(exp(- eta)))
    logq <- ifelse(eta < 0, - log1p(exp(eta)), - eta - log1p(exp(- eta)))
    logl <- sum(logp[y == 1]) + sum(logq[y == 0])
    return(logl - sum(beta^2) / 8)
}


###################################################
### code chunk number 4: metropolis-try-1
###################################################
set.seed(42)    # to get reproducible results
beta.init <- as.numeric(coefficients(out))
out <- metrop(lupost, beta.init, 1e3, x = x, y = y)
names(out)
out$accept


###################################################
### code chunk number 5: metropolis-try-2
###################################################
out <- metrop(out, scale = 0.1, x = x, y = y)
out$accept
out <- metrop(out, scale = 0.3, x = x, y = y)
out$accept
out <- metrop(out, scale = 0.5, x = x, y = y)
out$accept
out <- metrop(out, scale = 0.4, x = x, y = y)
out$accept


###################################################
### code chunk number 6: metropolis-try-3
###################################################
out <- metrop(out, nbatch = 1e4, x = x, y = y)
out$accept
out$time


###################################################
### code chunk number 7: fig1too
###################################################
plot(ts(out$batch))


###################################################
### code chunk number 8: fig1
###################################################
plot(ts(out$batch))


###################################################
### code chunk number 9: fig2too
###################################################
acf(out$batch)


###################################################
### code chunk number 10: fig2
###################################################
acf(out$batch)


###################################################
### code chunk number 11: metropolis-try-4
###################################################
out <- metrop(out, nbatch = 1e2, blen = 100,
    outfun = function(z, ...) c(z, z^2), x = x, y = y)
out$accept
out$time


###################################################
### code chunk number 12: metropolis-batch
###################################################
apply(out$batch, 2, mean)


###################################################
### code chunk number 13: metropolis-batch-too
###################################################
foo <- apply(out$batch, 2, mean)
mu <- foo[1:5]
sigmasq <- foo[6:10] - mu^2
mu
sigmasq


###################################################
### code chunk number 14: metropolis-mcse-mu
###################################################
mu.mcse <- apply(out$batch[ , 1:5], 2, sd) / sqrt(out$nbatch)
mu.mcse


###################################################
### code chunk number 15: metropolis-mcse-sigmasq
###################################################
u <- out$batch[ , 1:5]
v <- out$batch[ , 6:10]
ubar <- apply(u, 2, mean)
vbar <- apply(v, 2, mean)
deltau <- sweep(u, 2, ubar)
deltav <- sweep(v, 2, vbar)
foo <- sweep(deltau, 2, ubar, "*")
sigmasq.mcse <- sqrt(apply((deltav - 2 * foo)^2, 2, mean) / out$nbatch)
sigmasq.mcse


###################################################
### code chunk number 16: metropolis-mcse-sigmasq-too
###################################################
sqrt(mean(((v[ , 2] - vbar[2]) - 2 * ubar[2] * (u[ , 2] - ubar[2]))^2) /
    out$nbatch)


###################################################
### code chunk number 17: metropolis-mcse-sigma
###################################################
sigma <- sqrt(sigmasq)
sigma.mcse <- sigmasq.mcse / (2 * sigma)
sigma
sigma.mcse


###################################################
### code chunk number 18: metropolis-try-5
###################################################
out <- metrop(out, nbatch = 5e2, blen = 400, x = x, y = y)
out$accept
out$time
foo <- apply(out$batch, 2, mean)
mu <- foo[1:5]
sigmasq <- foo[6:10] - mu^2
mu
sigmasq
mu.mcse <- apply(out$batch[ , 1:5], 2, sd) / sqrt(out$nbatch)
mu.mcse
u <- out$batch[ , 1:5]
v <- out$batch[ , 6:10]
ubar <- apply(u, 2, mean)
vbar <- apply(v, 2, mean)
deltau <- sweep(u, 2, ubar)
deltav <- sweep(v, 2, vbar)
foo <- sweep(deltau, 2, ubar, "*")
sigmasq.mcse <- sqrt(apply((deltav - 2 * foo)^2, 2, mean) / out$nbatch)
sigmasq.mcse
sigma <- sqrt(sigmasq)
sigma.mcse <- sigmasq.mcse / (2 * sigma)
sigma
sigma.mcse


###################################################
### code chunk number 19: tab1
###################################################
foo <- rbind(mu, mu.mcse)
dimnames(foo) <- list(c("estimate", "MCSE"),
    c("constant", paste("$x_", 1:4, "$", sep = "")))
library(xtable)
print(xtable(foo, digits = rep(4, 6),
    align = c("l", rep("c", 5))), floating = FALSE,
    caption.placement = "top",
    sanitize.colnames.function = function(x) return(x))


###################################################
### code chunk number 20: tab1
###################################################
foo <- rbind(sigmasq, sigmasq.mcse)
dimnames(foo) <- list(c("estimate", "MCSE"),
    c("constant", paste("$x_", 1:4, "$", sep = "")))
library(xtable)
print(xtable(foo, digits = rep(4, 6),
    align = c("l", rep("c", 5))), floating = FALSE,
    caption.placement = "top",
    sanitize.colnames.function = function(x) return(x))


###################################################
### code chunk number 21: tab1
###################################################
foo <- rbind(sigma, sigma.mcse)
dimnames(foo) <- list(c("estimate", "MCSE"),
    c("constant", paste("$x_", 1:4, "$", sep = "")))
library(xtable)
print(xtable(foo, digits = rep(4, 6),
    align = c("l", rep("c", 5))), floating = FALSE,
    caption.placement = "top",
    sanitize.colnames.function = function(x) return(x))


###################################################
### code chunk number 22: time
###################################################
cat(out$time[1], "\n")


###################################################
### code chunk number 23: x
###################################################
n <- 2e4
rho <- 0.99
x <- arima.sim(model = list(ar = rho), n = n)


###################################################
### code chunk number 24: figgamtoo
###################################################
out <- initseq(x)
plot(seq(along = out$Gamma.pos) - 1, out$Gamma.pos,
        xlab = "k", ylab = expression(Gamma[k]), type = "l")
lines(seq(along = out$Gamma.dec) - 1, out$Gamma.dec, lty = "dotted")
lines(seq(along = out$Gamma.con) - 1, out$Gamma.con, lty = "dashed")


###################################################
### code chunk number 25: figgam
###################################################
out <- initseq(x)
plot(seq(along = out$Gamma.pos) - 1, out$Gamma.pos,
        xlab = "k", ylab = expression(Gamma[k]), type = "l")
lines(seq(along = out$Gamma.dec) - 1, out$Gamma.dec, lty = "dotted")
lines(seq(along = out$Gamma.con) - 1, out$Gamma.con, lty = "dashed")


###################################################
### code chunk number 26: assvar
###################################################
out$var.con
(1 + rho) / (1 - rho) * 1 / (1 - rho^2)


###################################################
### code chunk number 27: batx
###################################################
blen <- 5
x.batch <- apply(matrix(x, nrow = blen), 2, mean)
bout <- initseq(x.batch)


###################################################
### code chunk number 28: figgambattoo
###################################################
plot(seq(along = bout$Gamma.con) - 1, bout$Gamma.con,
        xlab = "k", ylab = expression(Gamma[k]), type = "l")


###################################################
### code chunk number 29: figgambat
###################################################
plot(seq(along = bout$Gamma.con) - 1, bout$Gamma.con,
        xlab = "k", ylab = expression(Gamma[k]), type = "l")


###################################################
### code chunk number 30: compvar
###################################################
out$var.con
bout$var.con * blen


###################################################
### code chunk number 31: ci-con
###################################################
mean(x) + c(-1, 1) * qnorm(0.975) * sqrt(out$var.con / length(x))
mean(x.batch) + c(-1, 1) * qnorm(0.975) * sqrt(bout$var.con / length(x.batch))


