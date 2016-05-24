#  Copyright (C) 2012 Yohan Chalabi
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License as
#  published by the Free Software Foundation; either version 2 or 3 of
#  the License (at your option).
#
#  This program is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

test.qdgl <- function() {

    Rqdgl <- function(p, lambdas) {

        lambda1 <- lambdas[1]
        lambda2 <- lambdas[2]
        lambda3 <- lambdas[3]
        lambda4 <- lambdas[4]

        stopifnot(lambda2 > 0)

        q <- numeric(length(p))

        for (i in seq_along(p)) {

            u <- p[i]

            ## x <-
            ##     if (u > 1 || u < 0) {
            ##         NaN
            ##     } else if (u == 0) {
            ##         if (lambda3 < 1)
            ##             Inf
            ##         else
            ##             1 / lambda2
            ##     } else if (u == 1) {
            ##         if (lambda4 < 1)
            ##             Inf
            ##         else
            ##             1 / lambda2
            ##     } else if (lambda3 == 0 && lambda4 == 0) {
            ##         1 / (u - u^2) / lambda2
            ##     } else if (u == 0 && lambda3 > 0) {
            ##         0
            ##     } else if (lambda4 == 0 && lambda3 != 0) {
            ##         (1 / (1 - u) + u^(lambda3 - 1)) / lambda2
            ##     } else if (lambda3 == 0 && lambda4 != 0) {
            ##         (u - 1 - (1 - u)^lambda4 * u) / ((u - 1) * u * lambda2)
            ##     } else {
            ##         (u^(lambda3 - 1) + (1 - u)^(lambda4 - 1)) / lambda2
            ##     }

            x <-
                if (u > 1 || u < 0) {
                    NaN
                } else {
                    (u^(lambda3 - 1) + (1 - u)^(lambda4 - 1)) / lambda2
                }

            q[i] <- x
        }

        q
    }

    tol <- 5 * .Machine$double.eps

    ###################################
    # generate different sets of parameter values
    iqr <- c(1, 2)
    med <- c(-1, 0, 1)
    chi <- c(-.75, -.5, -.25, 0, .25, .5, .75)    # take care to include
    xi <- unique(sort(.5 * (1. + c(chi, -chi))))  # special cases
    params <- expand.grid(med, iqr, chi, xi)

    params <- rbind(params,
                    c(0, 1, 0, .5 - 1/sqrt(5)), # special case
                    c(0, 1, 0, .5 - 2/sqrt(17)),
                    c(0, 1, 0, .5),
                    c(0, 1,  1, 0),
                    c(0, 1, -1, 0))

    # special case with different location and scale
    params <- rbind(params,
                    c(-1, 2, 0, .5 - 1/sqrt(5)),
                    c(-1, 2, 0, .5 - 2/sqrt(17)),
                    c(-1, 2, 0, .5),
                    c(-1, 2,  1, 0),
                    c(-1, 2, -1, 0))
    params <- rbind(params,
                    c(1, 2, 0, .5 - 1/sqrt(5)),
                    c(1, 2, 0, .5 - 2/sqrt(17)),
                    c(1, 2, 0, .5),
                    c(1, 2,  1, 0),
                    c(1, 2, -1, 0))

    ###################################
    # compare both implementation
    p <- seq(0, 1, len = 1e3)
    for (i in seq.int(nrow(params))) {
        pars <- as.numeric(params[i, ])
        lambdas <- CSW2FKML(pars[1], pars[2], pars[3], pars[4])
        if (!(i %in% c(295, 300, 305)))
        checkEquals(Rqdgl(p, lambdas), qdgl(p, pars), tol)
    }

    # deal with special case which correspond to the included uniform
    # distribuion. Rqdgl fails at p = 0,1 because of numerical errors
    # in the transformed lambdas.
    pars <- as.numeric(params[295, ])
    checkEquals(Rqdgl(p, c(0, 1, 1, 1)), qdgl(p, pars))
    pars <- as.numeric(params[300, ])
    checkEquals(Rqdgl(p, c(-1, .5, 1, 1)), qdgl(p, pars))
    pars <- as.numeric(params[305, ])
    checkEquals(Rqdgl(p, c(1, .5, 1, 1)), qdgl(p, pars))

    ###################################
    # test included distribution

    p <- seq(0, 1, len = 1e3)

    # uniform
    a <- 1; b <- 5
    pars <- c(.5 * (a + b), .5 * (b - a), 0, .5 - 1/sqrt(5))
    qd <- qdgl(p, pars)
    # -> note it fails for the boudary values because of numerical
    # -> impresition when setting the critical parameters
    checkEquals(rep(4., length(p)-2), qd[-c(1, length(p))], tol)
    pars <- c(.5 * (a + b), .5 * (b - a), 0, .5 - 2/sqrt(17))
    qd <- qdgl(p, pars)
    checkEquals(rep(4., length(p)), qd, tol)

    # logistic
    qdlogis <- function(p, loc, scale)
        ifelse (p == 0, Inf , scale / (1/p - 1) / p^2)
    loc <- 2; scale <- 3
    pars <- c(loc, scale * log(9), 0, .5)
    checkEquals(qdlogis(p, loc, scale), qdgl(p, pars), tol)

    # exponential
    qdexp <- function(p, rate)
        ifelse(p == 0, 1 / rate, ifelse(p == 1, Inf, 1 / (1 - p) / rate))
    rate <- 3

    pars <- c(log(2)/rate, log(3)/rate, 1, 0)
    checkEquals(qdexp(p, rate), qdgl(p, pars), tol)

    pars <- c(-log(2)/rate, log(3)/rate, -1, 0)
    checkEquals(qdexp(p, rate), qdgl(1-p, pars), tol)

    ###################################
    # test with NAs in vector of probabilities
    for (i in 1:nrow(params)) {
        pars <- as.numeric(params[i, ])
        p <- c(0, NA, .25, NA, .5, NA, .75, NA, 1, NA)
        q <- qdgl(p, pars)
        checkIdentical(q[c(FALSE, TRUE)], rep(NA_real_, 5))
    }

    ###################################
    # test with NaNs in vector of probabilities
    for (i in 1:nrow(params)) {
        pars <- as.numeric(params[i, ])
        p <- c(0, NaN, .25, NaN, .5, NaN, .75, NaN, 1, NaN)
        q <- qdgl(p, pars)
        checkIdentical(q[c(FALSE, TRUE)], rep(NaN, 5))
    }

    ###################################
    # test with values larger than one or smaller than 0

    oo <- options(warn = -1)
    on.exit(options(oo))

    for (i in 1:nrow(params)) {
        pars <- as.numeric(params[i, ])
        p <- c(-1, .5, 1.5)
        q <- qdgl(p, pars)
        checkIdentical(q[c(1, 3)], rep(NaN, 2))
        checkIdentical(q[2], qdgl(p[2], pars))
    }

    options(warn = 2)

    for (i in 1:nrow(params)) {
        pars <- as.numeric(params[i, ])
        p <- c(-1, .5, 1.5)
        checkException(qdgl(p, pars), silent = TRUE)
    }

    ###################################
    # test invalid options that should produce NaNs and warnings

    options(warn = -1)

    checkTrue(all(is.nan(qdgl(c(0, 0.5, 1), iqr = 0))))
    checkTrue(all(is.nan(qdgl(c(0, 0.5, 1), iqr = -1))))
    checkTrue(all(is.nan(qdgl(c(0, 0.5, 1), chi = -1))))
    checkTrue(all(is.nan(qdgl(c(0, 0.5, 1), chi = 1))))
    checkTrue(all(is.nan(qdgl(c(0, 0.5, 1), chi = 2))))
    checkTrue(all(is.nan(qdgl(c(0, 0.5, 1), chi = -2))))
    checkTrue(all(is.nan(qdgl(c(0, 0.5, 1), xi = 0))))
    checkTrue(all(is.nan(qdgl(c(0, 0.5, 1), xi = 1))))
    checkTrue(all(is.nan(qdgl(c(0, 0.5, 1), xi = -1))))
    checkTrue(all(is.nan(qdgl(c(0, 0.5, 1), xi = 2))))

    options(warn = 2)

    checkException(qdgl(c(0, 0.5, 1), iqr = 0), silent = TRUE)
    checkException(qdgl(c(0, 0.5, 1), iqr = -1), silent = TRUE)
    checkException(qdgl(c(0, 0.5, 1), chi = -1), silent = TRUE)
    checkException(qdgl(c(0, 0.5, 1), chi = 1), silent = TRUE)
    checkException(qdgl(c(0, 0.5, 1), chi = 2), silent = TRUE)
    checkException(qdgl(c(0, 0.5, 1), chi = -2), silent = TRUE)
    checkException(qdgl(c(0, 0.5, 1), xi = 0), silent = TRUE)
    checkException(qdgl(c(0, 0.5, 1), xi = 1), silent = TRUE)
    checkException(qdgl(c(0, 0.5, 1), xi = -1), silent = TRUE)
    checkException(qdgl(c(0, 0.5, 1), xi = 2), silent = TRUE)

}
