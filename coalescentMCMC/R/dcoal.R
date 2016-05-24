## dcoal.R (2013-12-03)

##   pdf of Various Time-Dependent Coalescent Models

## Copyright 2012-2013 Emmanuel Paradis

## This file is part of the R-package `coalescentMCMC'.
## See the file ../COPYING for licensing issues.

dcoal.step <- function(phy, theta0, theta1, tau, log = FALSE)
{
    ## if theta0 = theta1 use dcoal():
    if (theta0 == theta1) return(dcoal(phy, theta0, log))
    ## branching times from the most recent to the oldest:
    x <- c(0, sort(branching.times(phy)))
    n <- length(x)
    f <- F <- numeric(n)
    s <- x <= tau
    f[s] <- theta0 # theta(t) for x <= tau
    f[!s] <- theta1 # theta(t) for x > tau
    ## the primitive of '1/f' is 'x/theta0' for x <= tau...
    F[s] <- x[s]/theta0
    ## ... see paper for x > tau:
    F[!s] <- tau/theta0 + (x[!s] - tau)/theta1
    Ncomb <- choose(n:2, 2)
    p <- sum(log(Ncomb) - log(f[-1]) - 2 * Ncomb * (F[-1] - F[-n]))
    if (!log) p <- exp(p)
    p
}

dcoal.linear <- function(phy, theta0, thetaT, TMRCA, log = FALSE)
{
    if (theta0 == thetaT) return(dcoal(phy, theta0, log))
    ## branching times from the most recent to the oldest:
    x <- c(0, sort(branching.times(phy)))
    kappa <- (thetaT - theta0)/TMRCA
    f <- theta0 + x*kappa # = theta(t)
    ## don't need to compute the inverse of theta(t)
    n <- length(f)
    Ncomb <- choose(n:2, 2)

    ## the primitive of 1/f is log(f)/f'
    lnf <- log(f)
    p <- sum(log(Ncomb) - lnf[-1] - 2 * Ncomb * (lnf[-1] - lnf[-n])/kappa)
    if (!log) p <- exp(p)
    p
}

dcoal.time2 <- function(phy, theta0, rho1, rho2, tau, log = FALSE)
{
    ## if rho1 = rho2 = 0 use dcoal():
    if (rho1 == rho2) {
        if (!rho1) return(dcoal(phy, theta0, log))
        else return(dcoal.time(phy, theta0, rho1, log))
    }
    ## branching times from the most recent to the oldest:
    x <- c(0, sort(branching.times(phy)))
    ## \theta(t) = \theta_0 e^{\rho_1 t} \quad t \le \tau
    ## \theta(t) = \theta_0 e^{\rho_1 \tau} e^{\rho_2 t} \quad t \gt \tau
    ##           = \theta_0 e^{\rho_1 \tau + \rho_2 t}
    n <- length(x)
    f <- F <- numeric(n)
    s <- x <= tau
    f[s] <- exp(-rho1*x[s])/theta0 # inverse of theta(t) for x <= tau
    f[!s] <- exp(-rho2*x[!s] - (rho1 - rho2)*tau)/theta0 # inverse of theta(t) for x > tau
    ## the primitive of 'f' is '-(f - 1/theta0)/rho1' for x <= tau...
    F[s] <- -(f[s] - 1/theta0)/rho1
    ## ... see paper for x > tau:
    A <- exp(-rho1*tau)/theta0
    F[!s] <- -(A - 1/theta0)/rho1 - (f[!s] - A)/rho2
    Ncomb <- choose(n:2, 2)
    p <- sum(log(Ncomb) + log(f[-1]) - 2 * Ncomb * (F[-1] - F[-n]))
    if (!log) p <- exp(p)
    p
}

dcoal.time <- function(phy, theta0, rho, log = FALSE)
{
    ## returns NaN if rho = 0; use dcoal() instead:
    if (!rho) return(dcoal(phy, theta0, log))
    ## branching times from the most recent to the oldest:
    x <- c(0, sort(branching.times(phy)))
    ## \theta(t) = \theta_0 e^{\rho t}
    f <- exp(-rho * x)/theta0 # inverse of theta(t)
    n <- length(f)
    Ncomb <- choose(n:2, 2)
    ## the primitive of 'f' is '-f/rho'
    p <- sum(log(Ncomb) + log(f[-1]) + 2 * Ncomb * (f[-1] - f[-n])/rho)

    gr.theta0 <- -(n - 1)/theta0 - sum(2 * Ncomb * (f[-1] - f[-n])/(theta0 * rho))
    gr.rho <- sum(-x[-1] + 2 * Ncomb * (-(f[-1] - f[-n])/rho^2 + (x[-1]*f[-1] - x[-n]*f[-n])/rho))
    if (log) attr(p, "gradient") <- c(gr.theta0, gr.rho)
    else p <- exp(p)
    p
}

dcoal <- function(phy, theta, log = FALSE)
### this function is vectorized on 'theta'
{
    ## coalescent intervals from the oldest to most recent one:
    x <- rev(diff(c(0, sort(branching.times(phy)))))
    k <- 2:length(phy$tip.label)
    tmp <- k * (k - 1)/2 # choose(k, 2)
    tmp2 <- sum(x * tmp)
    sltmp <- sum(log(tmp))
    p <- sltmp - length(k) * log(theta) - 2 * tmp2/theta
    if (!log) p <- exp(p)
    p
}
