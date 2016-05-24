### R code from vignette source 'trunc.Rnw'

###################################################
### code chunk number 1: foo
###################################################
options(keep.source = TRUE, width = 70)


###################################################
### code chunk number 2: foo-pois
###################################################
k <- 2
theta <- seq(-100, 100, 10)
psi <- function(theta) {
    mu <- exp(theta)
    mu + ppois(k, lambda = mu, lower.tail = FALSE, log.p = TRUE)
}
tau <- function(theta) {
    mu <- exp(theta)
    beeta <- ppois(k + 1, lambda = mu, lower.tail = FALSE) /
        dpois(k + 1, lambda = mu)
    mu + (k + 1) / (beeta + 1)
}
qux <- function(theta) {
    mu <- exp(theta)
    beeta <- ppois(k + 1, lambda = mu, lower.tail = FALSE) /
        dpois(k + 1, lambda = mu)
    pbeeta <- ifelse(beeta < 1, beeta / (1 + beeta), 1 / (1 / beeta + 1))
    mu * (1 - (k + 1) / (beeta + 1) * (1 - (k + 1) / mu * pbeeta))
}


###################################################
### code chunk number 3: foo-pois-deriv
###################################################
epsilon <- 1e-6
pfoo <- tau(theta)
pbar <- (psi(theta + epsilon) - psi(theta)) / epsilon
all.equal(pfoo, pbar, tolerance = 10 * epsilon)
pfoo <- qux(theta)
pbar <- (tau(theta + epsilon) - tau(theta)) / epsilon
all.equal(pfoo, pbar, tolerance = 10 * epsilon)


###################################################
### code chunk number 4: foo-pois-mean
###################################################
theta <- seq(log(0.01), log(100), length = 51)
pfoo <- tau(theta)
pqux <- qux(theta)
pbar <- double(length(theta))
pbaz <- double(length(theta))
for (i in seq(along = theta)) {
    mu <- exp(theta[i])
    xxx <- seq(0, 10000)
    ppp <- dpois(xxx, lambda = mu)
    ppp[xxx <= k] <- 0
    ppp <- ppp / sum(ppp)
    pbar[i] <- sum(xxx * ppp)
    pbaz[i] <- sum((xxx - pbar[i])^2 * ppp)
}
all.equal(pfoo, pbar)
all.equal(pqux, pbaz)


###################################################
### code chunk number 5: foo-neg-bin
###################################################
k <- 2
alpha <- 2.22
mu <- 10^seq(-2, 2, 0.1)
theta <- log(mu) - log(mu + alpha)
psi <- function(theta) {
    stopifnot(all(theta < 0))
    mu <- (- alpha * exp(theta) / expm1(theta))
    alpha * log1p(1 / expm1(- theta)) +
        pnbinom(k, size = alpha, mu = mu, lower.tail = FALSE, log.p = TRUE)
}
tau <- function(theta) {
    stopifnot(all(theta < 0))
    mu <- (- alpha * exp(theta) / expm1(theta))
    p <- alpha / (mu + alpha)
    beetaup <- pnbinom(k + 1, size = alpha, mu = mu, lower.tail = FALSE)
    beetadn <- dnbinom(k + 1, size = alpha, mu = mu)
    beeta <- beetaup / beetadn
    beeta[beetaup == 0] <- 0
    result <- mu + (k + 1) / (beeta + 1) / p
    result[p == 0] <- Inf
    return(result)
}
qux <- function(theta) {
    stopifnot(all(theta < 0))
    mu <- (- alpha * exp(theta) / expm1(theta))
    p <- alpha / (mu + alpha)
    omp <- mu / (mu + alpha)
    beetaup <- pnbinom(k + 1, size = alpha, mu = mu, lower.tail = FALSE)
    beetadn <- dnbinom(k + 1, size = alpha, mu = mu)
    beeta <- beetaup / beetadn
    beeta[beetaup == 0] <- 0
    pbeeta <- ifelse(beeta < 1, beeta / (1 + beeta), 1 / (1 / beeta + 1))
    alpha * omp / p^2 - (k + 1) / p^2 / (1 + beeta) * ( - omp +
        (k + 1 + alpha) * omp / (1 + beeta) +
        ( alpha - p * (k + 1 + alpha) ) * pbeeta )
}


###################################################
### code chunk number 6: foo-neg-bin-deriv
###################################################
epsilon <- 1e-6
pfoo <- tau(theta)
pbar <- (psi(theta + epsilon) - psi(theta)) / epsilon
all.equal(pfoo, pbar, tolerance = 20 * epsilon)
pfoo <- qux(theta)
pbar <- (tau(theta + epsilon) - tau(theta)) / epsilon
all.equal(pfoo, pbar, tolerance = 40 * epsilon)


###################################################
### code chunk number 7: foo-neg-bin-mean
###################################################
pfoo <- tau(theta)
pqux <- qux(theta)
pbar <- double(length(theta))
pbaz <- double(length(theta))
for (i in seq(along = theta)) {
    mu <- (- alpha * exp(theta[i]) / expm1(theta[i]))
    xxx <- seq(0, 10000)
    ppp <- dnbinom(xxx, size = alpha, mu = mu)
    ppp[xxx <= k] <- 0
    ppp <- ppp / sum(ppp)
    pbar[i] <- sum(xxx * ppp)
    pbaz[i] <- sum((xxx - pbar[i])^2 * ppp)
}
all.equal(pfoo, pbar)
all.equal(pqux, pbaz)


###################################################
### code chunk number 8: valid-neg-bin
###################################################
alpha <- 2.22
p <- 0.5
k <- 20
m <- max(ceiling((k + 1 + alpha) * p - alpha), 0)
m

nsim <- 1e6
y <- rnbinom(nsim, size = alpha + m, prob = p)
xprop <- y + m
aprop <- exp(lfactorial(y) - lfactorial(xprop) + lfactorial(k + 1) -
    lfactorial(k + 1 - m)) * as.numeric(xprop > k)
max(aprop)
x <- xprop[runif(nsim) < aprop]
n <- length(x)

fred <- tabulate(x)
xfred <- seq(along = fred)
pfred <- dnbinom(xfred, size = alpha, prob = p)
pfred[xfred <= k] <- 0
pfred <- pfred / sum(pfred)
mfred <- max(xfred[n * pfred > 5])
o <- fred
o[mfred] <- sum(o[seq(mfred, length(fred))])
o <- o[seq(k + 1, mfred)]
e <- n * pfred
e[mfred] <- sum(e[seq(mfred, length(fred))])
e <- e[seq(k + 1, mfred)]
chisqstat <- sum((o - e)^2 / e)
pchisq(chisqstat, lower.tail = FALSE, df = length(o))


###################################################
### code chunk number 9: perf-neg-bin
###################################################
length(x) / nsim
rho <- function(m, p) {
    exp(lfactorial(k + 1) - lfactorial(k + 1 - m) +
    lgamma(alpha) - lgamma(m + alpha) + m * (log(p) - log1p(- p)) +
    pnbinom(k, size = alpha, prob = p, lower.tail = FALSE, log.p = TRUE))
}
rho(m, p)


###################################################
### code chunk number 10: perf-curve-neg-bin
###################################################
mu <- seq(0.01, k + 5, 0.01)
p <- alpha / (alpha + mu)
m <- pmax(ceiling((k + 1 + alpha) * p - alpha), 0)
r <- rho(m, p)


###################################################
### code chunk number 11: fig1
###################################################
plot(mu, r, xlab = expression(mu), ylab = "acceptance rate",
    type = "l")


###################################################
### code chunk number 12: valid-pois
###################################################
mu <- 2.22
k <- 20
m <- max(ceiling(k + 1 - mu), 0)
m

nsim <- 1e6
y <- rpois(nsim, lambda = mu)
xprop <- y + m
aprop <- exp(lfactorial(y) - lfactorial(xprop) + lfactorial(k + 1) -
    lfactorial(k + 1 - m)) * as.numeric(xprop > k)
max(aprop)
x <- xprop[runif(nsim) < aprop]
n <- length(x)

fred <- tabulate(x)
xfred <- seq(along = fred)
pfred <- dpois(xfred, lambda = mu)
pfred[xfred <= k] <- 0
pfred <- pfred / sum(pfred)
mfred <- max(xfred[n * pfred > 5])
o <- fred
o[mfred] <- sum(o[seq(mfred, length(fred))])
o <- o[seq(k + 1, mfred)]
e <- n * pfred
e[mfred] <- sum(e[seq(mfred, length(fred))])
e <- e[seq(k + 1, mfred)]
chisqstat <- sum((o - e)^2 / e)
pchisq(chisqstat, lower.tail = FALSE, df = length(o))


###################################################
### code chunk number 13: perf-pois
###################################################
length(x) / nsim
rho <- function(m, mu) {
    exp(lfactorial(k + 1) - lfactorial(k + 1 - m) - m * log(mu) +
    ppois(k, lambda = mu, lower.tail = FALSE, log.p = TRUE))
}
rho(m, mu)


###################################################
### code chunk number 14: perf-curve-pois
###################################################
mu <- seq(0.01, k + 5, 0.01)
m <- pmax(ceiling(k + 1 - mu), 0)
r <- rho(m, mu)


###################################################
### code chunk number 15: fig2
###################################################
plot(mu, r, xlab = expression(mu), ylab = "acceptance rate",
    type = "l")


###################################################
### code chunk number 16: perf-curve-pois-various
###################################################
kseq <- c(0, 1, 2, 20, 100)
mseq <- double(length(kseq))
for (i in seq(along = kseq)) {
    k <- kseq[i]
    mu <- seq(0.01, k + 5, 0.01)
    m <- pmax(ceiling(k + 1 - mu), 0)
    r <- rho(m, mu)
    mseq[i] <- min(r)
}


