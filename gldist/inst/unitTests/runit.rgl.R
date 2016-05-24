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

test.rgl <- function () {

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

    tol <- 5 * .Machine$double.eps

    ###################################
    # check that we have the same sequence when calling rgl or using
    # qgl and runif.
    for (i in seq.int(nrow(params))) {
        pars <- as.numeric(params[i, ])
        set.seed(2012)
        x1 <- rgl(1000L, pars)
        set.seed(2012)
        x2 <- qgl(runif(1000L), pars)
        checkEquals(x1, x2, tol = tol)
    }

    ###################################
    # test invalid options that should produce NaNs and warnings

    oo <- options(warn = -1)
    on.exit(options(oo))

    checkTrue(all(is.nan(rgl(3, iqr = 0))))
    checkTrue(all(is.nan(rgl(3, iqr = -1))))
    checkTrue(all(is.nan(rgl(3, chi = -1))))
    checkTrue(all(is.nan(rgl(3, chi = 1))))
    checkTrue(all(is.nan(rgl(3, chi = 2))))
    checkTrue(all(is.nan(rgl(3, chi = -2))))
    checkTrue(all(is.nan(rgl(3, xi = 0))))
    checkTrue(all(is.nan(rgl(3, xi = 1))))
    checkTrue(all(is.nan(rgl(3, xi = -1))))
    checkTrue(all(is.nan(rgl(3, xi = 2))))

    options(warn = 2)

    checkException(rgl(3, iqr = 0), silent = TRUE)
    checkException(rgl(3, iqr = -1), silent = TRUE)
    checkException(rgl(3, chi = -1), silent = TRUE)
    checkException(rgl(3, chi = 1), silent = TRUE)
    checkException(rgl(3, chi = 2), silent = TRUE)
    checkException(rgl(3, chi = -2), silent = TRUE)
    checkException(rgl(3, xi = 0), silent = TRUE)
    checkException(rgl(3, xi = 1), silent = TRUE)
    checkException(rgl(3, xi = -1), silent = TRUE)
    checkException(rgl(3, xi = 2), silent = TRUE)

    options(oo)

}
