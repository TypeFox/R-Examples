## SimInf, a framework for stochastic disease spread simulations
## Copyright (C) 2015 - 2016  Stefan Engblom
## Copyright (C) 2015 - 2016  Stefan Widgren
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

library(SimInf)

## For debugging
sessionInfo()

## Define a tolerance
tol = 1e-8

## Local model parameters
l <- matrix(c(rep(91, 10), rep(182, 10), rep(273, 10), rep(365, 10)),
            nrow  = 4,
            byrow = TRUE)
storage.mode(l) <- "double"

## Distance matrix
d <- new("dgCMatrix",
         i = c(1L, 2L, 0L, 2L, 3L, 0L, 1L, 3L, 4L, 1L, 2L, 4L, 5L,
               2L, 3L, 5L, 6L, 3L, 4L, 6L, 7L, 4L, 5L, 7L, 8L, 5L,
               6L, 8L, 9L, 6L, 7L, 9L, 7L, 8L),
         p = c(0L, 2L, 5L, 9L, 13L, 17L, 21L, 25L, 29L, 32L, 34L),
         Dim = c(10L, 10L),
         Dimnames = list(NULL, NULL),
         x = c(1.4142135623731, 2.82842712474619, 1.4142135623731,
               1.4142135623731, 2.82842712474619, 2.82842712474619,
               1.4142135623731, 1.4142135623731, 2.82842712474619,
               2.82842712474619, 1.4142135623731, 1.4142135623731,
               2.82842712474619, 2.82842712474619, 1.4142135623731,
               1.4142135623731, 2.82842712474619, 2.82842712474619,
               1.4142135623731, 1.4142135623731, 2.82842712474619,
               2.82842712474619, 1.4142135623731, 1.4142135623731,
               2.82842712474619, 2.82842712474619, 1.4142135623731,
               1.4142135623731, 2.82842712474619, 2.82842712474619,
               1.4142135623731, 1.4142135623731, 2.82842712474619,
               1.4142135623731),
         factors = list())

## Check 'distance_matrix' method
d_obs <- distance_matrix(1:10, 1:10, 3)
stopifnot(is(d_obs, "dgCMatrix"))
stopifnot(identical(d_obs@i, d@i))
stopifnot(identical(d_obs@p, d@p))
stopifnot(all(abs(d_obs@x - d@x) < tol))

res <- tools::assertError(distance_matrix(1:10, 1:10, 3, "min_dist"))
stopifnot(length(grep("Invalid 'min_dist' argument. Please provide 'min_dist' > 0.",
                      res[[1]]$message)) > 0)

res <- tools::assertError(distance_matrix(1:10, 1:10, 3, c(1, 2)))
stopifnot(length(grep("Invalid 'min_dist' argument. Please provide 'min_dist' > 0.",
                      res[[1]]$message)) > 0)

res <- tools::assertError(distance_matrix(1:10, 1:10, 3, -1))
stopifnot(length(grep("Invalid 'min_dist' argument. Please provide 'min_dist' > 0.",
                      res[[1]]$message)) > 0)

## Check 'data' argument to C function 'siminf_ldata_sp'
res <- tools::assertError(
    .Call("siminf_ldata_sp", NULL, d, PACKAGE = "SimInf"))
stopifnot(length(grep("Invalid 'data' argument",
                      res[[1]]$message)) > 0)

res <- tools::assertError(
    .Call("siminf_ldata_sp", d, d, PACKAGE = "SimInf"))
stopifnot(length(grep("Invalid 'data' argument",
                      res[[1]]$message)) > 0)

res <- tools::assertError(
    .Call("siminf_ldata_sp", 1:10, d, PACKAGE = "SimInf"))
stopifnot(length(grep("Invalid 'data' argument",
                      res[[1]]$message)) > 0)

## Check 'distance' argument to C function 'siminf_ldata_sp'
res <- tools::assertError(
    .Call("siminf_ldata_sp", l, NULL, PACKAGE = "SimInf"))
stopifnot(length(grep("Invalid 'distance' argument",
                      res[[1]]$message)) > 0)

res <- tools::assertError(
    .Call("siminf_ldata_sp", l, l, PACKAGE = "SimInf"))
stopifnot(length(grep("Invalid 'distance' argument",
                      res[[1]]$message)) > 0)

## Check non-equal number of nodes in 'distance' and 'data'
res <- tools::assertError(
    .Call("siminf_ldata_sp", l[, -1], d, PACKAGE = "SimInf"))
stopifnot(length(grep("The number of nodes in 'data' and 'distance' are not equal",
                      res[[1]]$message)) > 0)

## Check 'ldata'
ldata_exp <- structure(c(91, 182, 273, 365, 1, 0.499999999999996, 2, 0.125,
                         -1, 0, 0, 0, 0, 0, 91, 182, 273, 365, 0,
                         0.499999999999996, 2, 0.499999999999996, 3, 0.125,
                         -1, 0, 0, 0, 91, 182, 273, 365, 0, 0.125, 1,
                         0.499999999999996, 3, 0.499999999999996, 4, 0.125,
                         -1, 0, 91, 182, 273, 365, 1, 0.125, 2,
                         0.499999999999996, 4, 0.499999999999996, 5, 0.125,
                         -1, 0, 91, 182, 273, 365, 2, 0.125, 3,
                         0.499999999999996, 5, 0.499999999999996, 6, 0.125,
                         -1, 0, 91, 182, 273, 365, 3, 0.125, 4,
                         0.499999999999996, 6, 0.499999999999996, 7, 0.125,
                         -1, 0, 91, 182, 273, 365, 4, 0.125, 5,
                         0.499999999999996, 7, 0.499999999999996, 8, 0.125,
                         -1, 0, 91, 182, 273, 365, 5, 0.125, 6,
                         0.499999999999996, 8, 0.499999999999996, 9, 0.125,
                         -1, 0, 91, 182, 273, 365, 6, 0.125, 7,
                         0.499999999999996, 9, 0.499999999999996, -1, 0, 0,
                         0, 91, 182, 273, 365, 7, 0.125, 8, 0.499999999999996,
                         -1, 0, 0, 0, 0, 0), .Dim = c(14L, 10L))

ldata_obs <- .Call("siminf_ldata_sp", l, d, PACKAGE = "SimInf")
stopifnot(all(abs(ldata_obs - ldata_exp) < tol))

## Check identical coordinates
res <- tools::assertError(
    distance_matrix(x = c(1,10,1), y = c(1,10,1), cutoff = 20))
stopifnot(length(grep("Identical coordinates. Please provide a minimum distance.",
                      res[[1]]$message)) > 0)

d_exp <- new("dgCMatrix",
             i = c(1L, 2L, 0L, 2L, 0L, 1L),
             p = c(0L, 2L, 4L, 6L),
             Dim = c(3L, 3L),
             Dimnames = list(NULL, NULL),
             x = c(12.7279220613579, 2, 12.7279220613579,
                   12.7279220613579, 2, 12.7279220613579),
         factors = list())
d_obs <- distance_matrix(x = c(1,10,1), y = c(1,10,1), cutoff = 20, min_dist = 2)
stopifnot(is(d_obs, "dgCMatrix"))
stopifnot(identical(d_obs@i, d_exp@i))
stopifnot(identical(d_obs@p, d_exp@p))
stopifnot(all(abs(d_obs@x - d_exp@x) < tol))
