### R code from vignette source 'morph.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: foo
###################################################
options(keep.source = TRUE, width = 60)
foo <- packageDescription("mcmc")


###################################################
### code chunk number 2: morph.Rnw:98-100
###################################################
library(mcmc)
h2 <- morph(b=1)


###################################################
### code chunk number 3: morph.Rnw:105-107
###################################################
lud <- function(x) dt(x, df=3, log=TRUE)
lud.induced <- h2$lud(lud)


###################################################
### code chunk number 4: morph.Rnw:110-114
###################################################
curve(exp(Vectorize(lud.induced)(x)), from = -3, to = 3, lty = 2,
    xlab = "t", ylab = "density")
curve(exp(lud(x)), add = TRUE)
legend("topright", c("t density", "induced density"), lty=1:2)


###################################################
### code chunk number 5: morph.Rnw:123-129
###################################################
lud(1:4)
lud(1)
foo <- try(lud.induced(1:4))
class(foo)
cat(foo, "\n")
lud.induced(1)


###################################################
### code chunk number 6: set-seed
###################################################
set.seed(42)


###################################################
### code chunk number 7: morph.Rnw:146-147
###################################################
out <- morph.metrop(lud, 0, blen=100, nbatch=100, morph=morph(b=1))


###################################################
### code chunk number 8: morph.Rnw:153-155
###################################################
# adjust scale to find a roughly 20% acceptance rate
out$accept


###################################################
### code chunk number 9: morph.Rnw:161-163
###################################################
out <- morph.metrop(out, scale=4)
out$accept


###################################################
### code chunk number 10: fig0too
###################################################
acf(out$batch)


###################################################
### code chunk number 11: fig0
###################################################
acf(out$batch)


###################################################
### code chunk number 12: morph.Rnw:187-188
###################################################
t.test(out$batch)


###################################################
### code chunk number 13: morph.Rnw:191-193
###################################################
colMeans(out$batch)
apply(out$batch, 2, sd) / sqrt(out$nbatch)


###################################################
### code chunk number 14: unmorph-metrop-adjust
###################################################
out.unmorph <- metrop(lud, 0, blen=1000, nbatch=1)
out.unmorph$accept
out.unmorph <- metrop(out.unmorph, scale=4)
out.unmorph$accept
out.unmorph <- metrop(out.unmorph, scale=6)
out.unmorph$accept


###################################################
### code chunk number 15: unmorph-metrop-t-long-run
###################################################
lout <- suppressWarnings(try(load("morph1.rda"), silent = TRUE))
if (inherits(lout, "try-error")) {
    out.unmorph <- metrop(out.unmorph, blen = 1e5, nbatch = 1e3)
    save(out.unmorph, file = "morph1.rda")
} else {
    .Random.seed <- out.unmorph$final.seed
}
out.unmorph$accept


###################################################
### code chunk number 16: fig4too
###################################################
foo <- as.vector(out.unmorph$batch)
qqnorm(foo)
qqline(foo)


###################################################
### code chunk number 17: fig4
###################################################
foo <- as.vector(out.unmorph$batch)
qqnorm(foo)
qqline(foo)


###################################################
### code chunk number 18: shapiro-wilk
###################################################
shapiro.test(foo)


###################################################
### code chunk number 19: morph-metrop-t-long-run
###################################################
lout <- suppressWarnings(try(load("morph2.rda"), silent = TRUE))
if (inherits(lout, "try-error")) {
    out.morph <- morph.metrop(out, blen = 1e5, nbatch = 1e3)
    save(out.morph, file = "morph2.rda")
} else {
    .Random.seed <- out.morph$final.seed
}
out.morph$accept


###################################################
### code chunk number 20: fig5too
###################################################
foo <- as.vector(out.morph$batch)
qqnorm(foo)
qqline(foo)


###################################################
### code chunk number 21: fig5
###################################################
foo <- as.vector(out.morph$batch)
qqnorm(foo)
qqline(foo)


###################################################
### code chunk number 22: shapiro-wilk
###################################################
shapiro.test(foo)


###################################################
### code chunk number 23: def-posterior-binom
###################################################
lud.binom <- function(beta, M, x, n) {
  MB <- M %*% beta
  sum(x * MB) - sum(n * log(1 + exp(MB)))
}


###################################################
### code chunk number 24: convert
###################################################
dat <- as.data.frame(UCBAdmissions)
dat.split <- split(dat, dat$Admit)
dat.split <- lapply(dat.split,
                    function(d) {
                      val <- as.character(d$Admit[1])
                      d["Admit"] <- NULL
                      names(d)[names(d) == "Freq"] <- val
                      d
                    })
dat <- merge(dat.split[[1]], dat.split[[2]])


###################################################
### code chunk number 25: build-model-matrix
###################################################
formula <- cbind(Admitted, Rejected) ~ (Gender + Dept)^2
mf <- model.frame(formula, dat)
M <- model.matrix(formula, mf)


###################################################
### code chunk number 26: morph.Rnw:396-398
###################################################
xi <- 0.30
nu <- 5


###################################################
### code chunk number 27: lud-binom
###################################################
lud.berkeley <- function(B)
  lud.binom(B, M, dat$Admitted + xi * nu, dat$Admitted + dat$Rejected + nu)


###################################################
### code chunk number 28: morph.Rnw:410-419
###################################################
berkeley.out <- morph.metrop(lud.berkeley, rep(0, ncol(M)), blen=1000,
                             nbatch=1, scale=0.1, morph=morph(p=3))
berkeley.out$accept
berkeley.out <- morph.metrop(berkeley.out, scale=0.05)
berkeley.out$accept
berkeley.out <- morph.metrop(berkeley.out, scale=0.02)
berkeley.out$accept
berkeley.out <- morph.metrop(berkeley.out, blen=10000)
berkeley.out$accept


###################################################
### code chunk number 29: morph.Rnw:422-423
###################################################
berkeley.out <- morph.metrop(berkeley.out, blen=1, nbatch=100000)


###################################################
### code chunk number 30: morph.Rnw:428-433
###################################################
beta <- setNames(colMeans(berkeley.out$batch), colnames(M))
MB <- M %*% beta
dat$p <- dat$Admitted / (dat$Admitted + dat$Rejected)
dat$p.post <- exp(MB) / (1 + exp(MB))
dat


###################################################
### code chunk number 31: calculate-posterior-probabilities
###################################################
posterior.probabilities <-
  t(apply(berkeley.out$batch, 1,
          function(r) {
            eMB <- exp(M %*% r)
            eMB / (1 + eMB)
          }))
quants <- apply(posterior.probabilities, 2, quantile, prob=c(0.05, 0.95))
quants.str <- matrix(apply(quants, 2,
                           function(r) sprintf("[%0.2f, %0.2f]", r[1], r[2])),
                     nrow=2, byrow=TRUE)



###################################################
### code chunk number 32: fig1
###################################################
x <- (0:5) * 2 + 1
plot(x[c(1, 6)] + 0.5 * c(-1, 1), 0:1,
     xlab="Department", ylab="Probability", xaxt="n", type="n")
axis(1, x, LETTERS[1:6])
for(i in 1:6) {
  lines((x[i]-0.25)*c(1, 1), quants[1:2, i], lwd=2, col="gray")
  lines((x[i] + 0.25) * c(1, 1), quants[1:2, i + 6], lwd=2, col="gray")
  points(x[i] + 0.25 * c(-1, 1), dat$p.post[i + c(0, 6)], pch=c("F", "M"))
}


###################################################
### code chunk number 33: cauchy-data
###################################################
n <- 15
mu0 <- 50
sigma0 <- 10
x <- rcauchy(n, mu0, sigma0)
round(sort(x), 1)


###################################################
### code chunk number 34: cauchy-log-unnormalized-posterior
###################################################
lup <- function(theta) {
    if (any(is.na(theta)))
        stop("NA or NaN in input to log unnormalized density function")
    mu <- theta[1]
    sigma <- theta[2]
    if (sigma <= 0) return(-Inf)
    if (any(! is.finite(theta))) return(-Inf)
    result <- sum(dcauchy(x, mu, sigma, log = TRUE)) - log(sigma)
    if (! is.finite(result)) {
        warning(paste("Oops!  mu = ", mu, "and sigma =", sigma))
    }
    return(result)
}


###################################################
### code chunk number 35: cauchy-robust
###################################################
mu.twiddle <- median(x)
sigma.twiddle <- IQR(x)
c(mu.twiddle, sigma.twiddle)


###################################################
### code chunk number 36: cauchy-posterior-mode
###################################################
oout <- optim(c(mu.twiddle, sigma.twiddle), lup,
    control = list(fnscale = -1), hessian = TRUE)
stopifnot(oout$convergence == 0)
mu.hat <- oout$par[1]
sigma.hat <- oout$par[2]
c(mu.hat, sigma.hat)


###################################################
### code chunk number 37: cauchy-hessian
###################################################
oout$hessian


###################################################
### code chunk number 38: cauchy-se
###################################################
sqrt(- 1 / diag(oout$hessian))


###################################################
### code chunk number 39: cauchy-doit
###################################################
moo <- morph(b = 0.5, r = 7, center = c(mu.hat, sigma.hat))
mout <- morph.metrop(lup, c(mu.hat, sigma.hat), 1e4,
    scale = 3, morph = moo)
mout$accept
mout <- morph.metrop(mout)


###################################################
### code chunk number 40: cfig1too
###################################################
acf(mout$batch)


###################################################
### code chunk number 41: cfig1
###################################################
acf(mout$batch)


###################################################
### code chunk number 42: cfig2too
###################################################
mu <- mout$batch[ , 1]
i <- seq(1, mout$nbatch, by = 15)
out.sub <- density(mu[i])
out <- density(mu, bw = out.sub$bw)
plot(out)


###################################################
### code chunk number 43: cfig2
###################################################
mu <- mout$batch[ , 1]
i <- seq(1, mout$nbatch, by = 15)
out.sub <- density(mu[i])
out <- density(mu, bw = out.sub$bw)
plot(out)


###################################################
### code chunk number 44: cfig3
###################################################
sigma <- mout$batch[ , 2]
out.sub <- density(sigma[i])
out <- density(sigma, bw = out.sub$bw)
plot(out)


