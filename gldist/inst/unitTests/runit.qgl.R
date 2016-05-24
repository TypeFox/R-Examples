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

test.qgl <- function () {

    Rqgl <- function(p, lambdas) {

        lambda1 <- lambdas[1]
        lambda2 <- lambdas[2]
        lambda3 <- lambdas[3]
        lambda4 <- lambdas[4]

        stopifnot(lambda2 > 0)

        q <- numeric(length(p))

        for (i in seq_along(p)) {

            u <- p[i]

            x <-
                if (0 < u && u < 1) {
                    if (lambda3 != 0 && lambda4 != 0) {
                        lambda1 +
                            (u^lambda3 - 1)/lambda3/lambda2 -
                                ((1 - u)^lambda4 - 1)/lambda4/lambda2
                    } else if (lambda3 == 0 && lambda4 != 0) {
                        lambda1 +
                            log(u)/lambda2 -
                                ((1 - u)^lambda4 - 1)/lambda4 / lambda2
                    } else if (lambda3 != 0 && lambda4 == 0) {
                        lambda1 +
                            (u^lambda3 - 1)/lambda3/lambda2 -
                                log(1 - u)/lambda2
                    } else {
                        lambda1 +
                            log(u)/lambda2 -
                                log(1 - u)/lambda2
                    }
                } else if (u == 0) {
                    if (lambda3 <= 0)
                        -Inf
                    else
                        lambda1 - 1/lambda3/lambda2
                } else if (u == 1) {
                    if (lambda4 <= 0)
                        Inf
                    else
                        lambda1 + 1/lambda4/lambda2
                } else {
                    NaN
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
        checkEquals(Rqgl(p, lambdas), qgl(p, pars), tol)
    }

    ###################################
    # test included distribution

    p <- seq(0, 1, len = 1e3)

    # uniform
    a <- 1; b <- 5
    q <- qunif(p, min = a, max = b)
    med <- .5 * (a + b); iqr <- .5 * (b - a); chi <- 0; xi <- .5 - 1/sqrt(5)
    checkEquals(q, qgl(p, med, iqr, chi, xi), tol)
    med <- .5 * (a + b); iqr <- .5 * (b - a); chi <- 0; xi <- .5 - 2/sqrt(17)
    checkEquals(q, qgl(p, med, iqr, chi, xi), tol)

    # logistic
    loc <- 2; scale <- 3
    q <- qlogis(p, loc, scale)
    med <- loc; iqr <- scale * log(9); chi <- 0; xi <- .5
    checkEquals(q, qgl(p, med, iqr, chi, xi), tol)

    # exponential
    rate <- 3
    q <- qexp(p, rate)
    med <- log(2)/rate; iqr <- log(3)/rate; chi <- 1; xi <- 0
    checkEquals(q, qgl(p, med, iqr, chi, xi), tol)
    med <- -log(2)/rate; iqr <- log(3)/rate; chi <- -1; xi <- 0
    checkEquals(q, -qgl(1-p, med, iqr, chi, xi), tol)

    ###################################
    # test with NAs in vector of probabilities
    for (i in 1:nrow(params)) {
        pars <- as.numeric(params[i, ])
        p <- c(0, NA, .25, NA, .5, NA, .75, NA, 1, NA)
        q <- qgl(p, pars)
        checkIdentical(q[c(FALSE, TRUE)], rep(NA_real_, 5))
    }

    ###################################
    # test with NaNs in vector of probabilities
    for (i in 1:nrow(params)) {
        pars <- as.numeric(params[i, ])
        p <- c(0, NaN, .25, NaN, .5, NaN, .75, NaN, 1, NaN)
        q <- qgl(p, pars)
        checkIdentical(q[c(FALSE, TRUE)], rep(NaN, 5))
    }

    ###################################
    # test with values larger than one or smaller than 0

    oo <- options(warn = -1)
    on.exit(options(oo))

    for (i in 1:nrow(params)) {
        pars <- as.numeric(params[i, ])
        p <- c(-1, .5, 1.5)
        q <- qgl(p, pars)
        checkIdentical(q[c(1, 3)], rep(NaN, 2))
        checkIdentical(q[2], qgl(p[2], pars))
    }

    options(warn = 2)

    for (i in 1:nrow(params)) {
        pars <- as.numeric(params[i, ])
        p <- c(-1, .5, 1.5)
        checkException(qgl(p, pars), silent = TRUE)
    }

    ###################################
    # test invalid options that should produce NaNs and warnings

    options(warn = -1)

    checkTrue(all(is.nan(qgl(c(0, 0.5, 1), iqr = 0))))
    checkTrue(all(is.nan(qgl(c(0, 0.5, 1), iqr = -1))))
    checkTrue(all(is.nan(qgl(c(0, 0.5, 1), chi = -1))))
    checkTrue(all(is.nan(qgl(c(0, 0.5, 1), chi = 1))))
    checkTrue(all(is.nan(qgl(c(0, 0.5, 1), chi = 2))))
    checkTrue(all(is.nan(qgl(c(0, 0.5, 1), chi = -2))))
    checkTrue(all(is.nan(qgl(c(0, 0.5, 1), xi = 0))))
    checkTrue(all(is.nan(qgl(c(0, 0.5, 1), xi = 1))))
    checkTrue(all(is.nan(qgl(c(0, 0.5, 1), xi = -1))))
    checkTrue(all(is.nan(qgl(c(0, 0.5, 1), xi = 2))))

    options(warn = 2)

    checkException(qgl(c(0, 0.5, 1), iqr = 0), silent = TRUE)
    checkException(qgl(c(0, 0.5, 1), iqr = -1), silent = TRUE)
    checkException(qgl(c(0, 0.5, 1), chi = -1), silent = TRUE)
    checkException(qgl(c(0, 0.5, 1), chi = 1), silent = TRUE)
    checkException(qgl(c(0, 0.5, 1), chi = 2), silent = TRUE)
    checkException(qgl(c(0, 0.5, 1), chi = -2), silent = TRUE)
    checkException(qgl(c(0, 0.5, 1), xi = 0), silent = TRUE)
    checkException(qgl(c(0, 0.5, 1), xi = 1), silent = TRUE)
    checkException(qgl(c(0, 0.5, 1), xi = -1), silent = TRUE)
    checkException(qgl(c(0, 0.5, 1), xi = 2), silent = TRUE)

    options(oo)

}
