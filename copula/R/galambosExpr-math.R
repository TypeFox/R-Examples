## Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.


`galambosCopula.expr` <-
structure(expression(cdf = (u1 * u2)^(1 - ((1 - log(u1)/log(u1 *
    u2))^(-alpha) + (log(u1)/log(u1 * u2))^(-alpha))^(-1/alpha)),
    pdf = (((log(u1) * (-log(u1) + log(u1 * u2)))/log(u1 * u2)^2)^alpha *
        (2 * log(u1) * ((1 - log(u1)/log(u1 * u2))^(-alpha) +
            (log(u1)/log(u1 * u2))^(-alpha))^(2/alpha) * (log(u1) -
            log(u1 * u2)) - log(u1 * u2)^2 + ((1 - log(u1)/log(u1 *
            u2))^(-alpha) + (log(u1)/log(u1 * u2))^(-alpha))^(1/alpha) *
            log(u1 * u2) * (1 + alpha + log(u1 * u2))) + ((1 -
        log(u1)/log(u1 * u2))^(-alpha) + (log(u1)/log(u1 * u2))^(-alpha))^(1/alpha) *
        ((1 - log(u1)/log(u1 * u2))^(2 * alpha) * (log(u1) -
            log(u1 * u2)) * (log(u1) * ((1 - log(u1)/log(u1 *
            u2))^(-alpha) + (log(u1)/log(u1 * u2))^(-alpha))^(1/alpha) -
            log(u1 * u2)) + log(u1) * (log(u1)/log(u1 * u2))^(2 *
            alpha) * (((1 - log(u1)/log(u1 * u2))^(-alpha) +
            (log(u1)/log(u1 * u2))^(-alpha))^(1/alpha) * (log(u1) -
            log(u1 * u2)) + log(u1 * u2))))/((u1 * u2)^((1 -
        log(u1)/log(u1 * u2))^(-alpha) + (log(u1)/log(u1 * u2))^(-alpha))^(-1/alpha) *
        log(u1) * ((1 - log(u1)/log(u1 * u2))^(-alpha) + (log(u1)/log(u1 *
        u2))^(-alpha))^(2/alpha) * ((1 - log(u1)/log(u1 * u2))^alpha +
        (log(u1)/log(u1 * u2))^alpha)^2 * (log(u1) - log(u1 *
        u2))), deriv1cdf = (u1 * u2)^(1 - ((1 - log(u1)/log(u1 *
        u2))^(-alpha) + (log(u1)/log(u1 * u2))^(-alpha))^(-1/alpha)) *
        ((1 - ((1 - log(u1)/log(u1 * u2))^(-alpha) + (log(u1)/log(u1 *
            u2))^(-alpha))^(-1/alpha))/u1 + ((-(alpha * (log(u1)/(u1 *
            log(u1 * u2)^2) - 1/(u1 * log(u1 * u2))) * (1 - log(u1)/log(u1 *
            u2))^(-1 - alpha)) - alpha * (-(log(u1)/(u1 * log(u1 *
            u2)^2)) + 1/(u1 * log(u1 * u2))) * (log(u1)/log(u1 *
            u2))^(-1 - alpha)) * ((1 - log(u1)/log(u1 * u2))^(-alpha) +
            (log(u1)/log(u1 * u2))^(-alpha))^(-1 - 1/alpha) *
            log(u1 * u2))/alpha)), .Names = c("cdf", "pdf", "deriv1cdf"
))
`galambosCopula.algr` <-
structure(expression(cdf = {
    .expr1 <- u1 * u2
    .expr4 <- log(u1)/log(.expr1)
    .expr6 <- -alpha
    .value <- .expr1^(1 - ((1 - .expr4)^.expr6 + .expr4^.expr6)^(-1/alpha))
    .grad <- array(0, c(length(.value), 1), list(NULL, c("s")))
    .grad[, "s"] <- 0
    attr(.value, "gradient") <- .grad
    .value
}, pdf = {
    .expr1 <- log(u1)
    .expr3 <- u1 * u2
    .expr4 <- log(.expr3)
    .expr7 <- .expr4^2
    .expr11 <- .expr1/.expr4
    .expr12 <- 1 - .expr11
    .expr13 <- -alpha
    .expr16 <- .expr12^.expr13 + .expr11^.expr13
    .expr18 <- .expr16^(2/alpha)
    .expr20 <- .expr1 - .expr4
    .expr24 <- .expr16^(1/alpha)
    .expr31 <- 2 * alpha
    .value <- ((.expr1 * (-.expr1 + .expr4)/.expr7)^alpha * (2 *
        .expr1 * .expr18 * .expr20 - .expr7 + .expr24 * .expr4 *
        (1 + alpha + .expr4)) + .expr24 * (.expr12^.expr31 *
        .expr20 * (.expr1 * .expr24 - .expr4) + .expr1 * .expr11^.expr31 *
        (.expr24 * .expr20 + .expr4)))/(.expr3^.expr16^(-1/alpha) *
        .expr1 * .expr18 * (.expr12^alpha + .expr11^alpha)^2 *
        .expr20)
    .grad <- array(0, c(length(.value), 1), list(NULL, c("s")))
    .grad[, "s"] <- 0
    attr(.value, "gradient") <- .grad
    .value
}, deriv1cdf = {
    .expr1 <- u1 * u2
    .expr2 <- log(u1)
    .expr3 <- log(.expr1)
    .expr4 <- .expr2/.expr3
    .expr5 <- 1 - .expr4
    .expr6 <- -alpha
    .expr9 <- .expr5^.expr6 + .expr4^.expr6
    .expr10 <- -1
    .expr13 <- 1 - .expr9^(.expr10/alpha)
    .expr18 <- .expr2/(u1 * .expr3^2)
    .expr20 <- 1/(u1 * .expr3)
    .expr23 <- .expr10 - alpha
    .value <- .expr1^.expr13 * (.expr13/u1 + (-(alpha * (.expr18 -
        .expr20) * .expr5^.expr23) - alpha * (-.expr18 + .expr20) *
        .expr4^.expr23) * .expr9^(.expr10 - 1/alpha) * .expr3/alpha)
    .grad <- array(0, c(length(.value), 1), list(NULL, c("s")))
    .grad[, "s"] <- 0
    attr(.value, "gradient") <- .grad
    .value
}), .Names = c("cdf", "pdf", "deriv1cdf"))
