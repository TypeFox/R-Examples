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

# convert the paramters from the different parametrization

.S_FKML <- function(p, lambda3, lambda4) {

    S <- numeric(length(p))

    for (i in seq_along(p)) {
        u <- p[i]

        S[i] <-
            if (0 < u && u < 1) {

                if (lambda3 != 0 && lambda4 != 0) {

                    (u^lambda3 - 1) / lambda3 - ((1 - u)^lambda4 - 1)/lambda4

                } else if ((lambda3 == 0 && lambda4 != 0)) {

                    log(u) - ((1- u)^lambda4 - 1)/lambda4

                } else if ((lambda3 != 0 && lambda4 == 0)) {

                    (u^lambda3 - 1) / lambda3 - log(1 - u)

                } else {

                    log(u) - log(1 - u)
                }

            } else if (u == 0) {

                if (lambda3 <= 0) -Inf else -1/lambda3

            } else if (u == 1) {

                if (lambda4 <= 0) -Inf else 1/lambda4

            } else {

                NaN

            }
    }

    # Result
    S
}

FKML2CSW <- function(lambda1, lambda2, lambda3, lambda4) {

    if (length(lambda1) == 4L) {
        par <- lambda1
        lambda1 <- par[1]
        lambda2 <- par[2]
        lambda3 <- par[3]
        lambda4 <- par[4]
    }

    med <- lambda1 + .S_FKML(1/2, lambda3, lambda4) / lambda2
    iqr <- (.S_FKML(3/4, lambda3, lambda4) - .S_FKML(1/4, lambda3, lambda4)) / lambda2

    chi <- (lambda3 - lambda4) / sqrt(1 + (lambda3 - lambda4)^2)
    xi <- .5 - (lambda3 + lambda4) / (2 * sqrt(1 + (lambda3 + lambda4)^2))

    c(med = med, iqr = iqr, chi = chi, xi = xi)
}

CSW2FKML <- function(med, iqr, chi, xi) {

    if (length(med) == 4L) {
        par <- med
        med <- par[1]
        iqr <- par[2]
        chi <- par[3]
        xi <- par[4]
    }

    if (iqr <= 0) {

        c(lambda1 = NaN, lambda2 = NaN, lambda3 = NaN, lambda4 = NaN)

    } else if (chi == 1 && xi == 0) {

        lambda3 <- Inf
        lambda4 <- 0
        lambda2 <- (.S_FKML(3/4, lambda3, lambda4) - .S_FKML(1/4, lambda3, lambda4)) / iqr
        lambda1 <- med - .S_FKML(1/2, lambda3, lambda4) / lambda2
        c(lambda1 = lambda1, lambda2 = lambda2, lambda3 = lambda3, lambda4 = lambda4)

    } else if (chi == -1 && xi == 0) {

        lambda3 <- 0
        lambda4 <- Inf
        lambda2 <- (.S_FKML(3/4, lambda3, lambda4) - .S_FKML(1/4, lambda3, lambda4)) / iqr
        lambda1 <- med - .S_FKML(1/2, lambda3, lambda4) / lambda2
        c(lambda1 = lambda1, lambda2 = lambda2, lambda3 = lambda3, lambda4 = lambda4)

    } else if (chi <= -1 || chi >= 1 || xi <= 0 || xi >= 1) {

        c(lambda1 = NaN, lambda2 = NaN, lambda3 = NaN, lambda4 = NaN)

    } else {

        alpha <- .5 * (.5 - xi) / sqrt(xi * (1 - xi))
        beta <- .5 * chi / sqrt(1 - chi^2)
        lambda3 <- alpha + beta
        lambda4 <- alpha - beta
        lambda2 <- (.S_FKML(3/4, lambda3, lambda4) - .S_FKML(1/4, lambda3, lambda4)) / iqr
        lambda1 <- med - .S_FKML(1/2, lambda3, lambda4) / lambda2
        c(lambda1 = lambda1, lambda2 = lambda2, lambda3 = lambda3, lambda4 = lambda4)

    }

}
