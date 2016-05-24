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

fAckley <- function (x) {
    n <- length(x)
    exp1 <-  exp(-0.2 * sqrt(1/n * sum(x^2)))
    exp2 <- exp(1/n * sum(cos(2 * pi * x)))
    -20 * exp1 - exp2 + 20 + exp(1)
}


fGriewank <- function (x) {
    s <- sum(x^2) / 4000
    p <- prod(cos(x / sqrt(seq_len(length(x)))))
    1 + s - p
}


fRastrigin <- function (x) {
    10 * length(x) + sum(x^2 - 10 * cos(2 * pi * x))
}


fRosenbrock <- function (x) {
    e <- function (i) 100 * (x[i+1] - x[i]^2)^2 + (1 - x[i])^2
    sum(sapply(seq_len(length(x) - 1), e))
}


fSphere <- function (x) {
    sum(x^2)
}


fSummationCancellation <- function (x) {
    s <- function (i) abs(sum(x[seq_len(i)]))
    -1 / (10^-5 + sum(sapply(seq_len(length(x)), s)))
}
