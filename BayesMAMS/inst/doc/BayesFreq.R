## ------------------------------------------------------------------------
library("BayesMAMS")
ssbayes(k=2, nu=1, q0=c(0, 0, 0), eta=0.95, zeta=0.90, deltastar=0.5, prec="known",
        crit="1")

## ------------------------------------------------------------------------
k <- 2
alloc <- sqrt(k)
nu <- 1
deltastar <- 0.5
alpha <- 0.05
power <- 0.90

## ------------------------------------------------------------------------
ssfreq_bon <- ((qnorm(1 - alpha/k) + qnorm(power)) / (sqrt(nu) * deltastar))^2 *
              (1 + 1/sqrt(k))
ceiling(c(sqrt(k) * ssfreq_bon, rep(ssfreq_bon, k)))

## ------------------------------------------------------------------------
library("mvtnorm")
rho <- 1 / (1 + alloc)
corr <- matrix(rho, k, k) + diag(1 - rho, k)
quan <- qmvnorm(0.95, mean=rep(0, k), corr=corr)$quantile
ssfreq_dun <- ((quan + qnorm(power)) / (sqrt(nu) * deltastar))^2 * (1 + 1/alloc)
ceiling(c(sqrt(k) * ssfreq_dun, rep(ssfreq_dun, k)))

## ----message=FALSE-------------------------------------------------------
library("MAMS")
pstar <- pnorm(deltastar / sqrt(2 * 1/nu))
mams(K=k, J=1, r=1, r0=alloc, p=pstar, p0=0.5)

## ------------------------------------------------------------------------
ssbayes(k=2, nu=1, q0=c(0, 0, 0), eta=0.95, zeta=0.90, deltastar=0.5, prec="known",
        crit="2")

## ------------------------------------------------------------------------
ssfreq_unadj <- ((qnorm(1 - alpha) + qnorm(power)) / (sqrt(nu) * deltastar))^2 *
                (1 + 1/sqrt(k))
ceiling(c(sqrt(k) * ssfreq_unadj, rep(ssfreq_unadj, k)))

