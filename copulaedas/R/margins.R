# copulaedas: Estimation of Distribution Algorithms Based on Copulas
# Copyright (C) 2011-2015 Yasser Gonzalez Fernandez
# Copyright (C) 2011-2015 Marta Soto Ortiz
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

# Normal (pnorm and qnorm defined in the stats package).

fnorm <- function (x, lower, upper) {
    list(mean = mean(x), sd = sd(x))
}


# Truncated normal (ptruncnorm and qtruncnorm defined in the truncnorm package).

ftruncnorm <- function (x, lower, upper) {
    list(a = lower, b = upper, mean = mean(x), sd = sd(x))
}


# Kernel-smoothed empirical margins.

fkernel <- function (x, lower, upper) {
    list(X = x, h = bw.nrd0(x))
}

pkernel <- function (q, X, h) {
    rank(q) / (length(q) + 1)
}

qkernel <- function (p, X, h) {
    eps <- .Machine$double.eps^0.5
    maxIter <- 100
    quantiles <- quantile(X, p, names = FALSE)
    n <- length(X)
    f <- function (x) sum(dnorm((x - X) / h)) / (n * h)
    F <- function (x) sum(pnorm((x - X) / h)) / n
    sapply(seq(along = p),
        function (i) {
            iter <- 0
            x <- quantiles[i]
            Fx <- F(x) - p[i]
            while (abs(Fx) > eps && iter <= maxIter) {
                fx <- f(x)
                if (is.finite(fx) && abs(fx) > eps) {
                    x <- x - Fx / fx
                    Fx <- F(x) - p[i]
                    iter <- iter + 1
                } else {
                    break
                }
            }
            x
        })
}


# Truncated kernel-smoothed empirical margins.

ftrunckernel <- function (x, lower, upper) {
    c(a = lower, b = upper, fkernel(x))
}

ptrunckernel <- function (q, a, b, X, h) {
    rank(q) / (length(q) + 1)
}

qtrunckernel <- function (p, a, b, X, h) {
    F <- function (x) sum(pnorm((x - X) / h)) / length(X)
    Q <- qkernel
    Fa <- F(a)
    Fb <- F(b)
    r <- Q(Fa + p*(Fb - Fa), X, h)
    pmin(pmax(r, a), b)
}
