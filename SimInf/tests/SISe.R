## SimInf, a framework for stochastic disease spread simulations
## Copyright (C) 2015  Pavol Bauer
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

## Check invalid u0
res <- tools::assertError(SISe(u0 = "u0"))
stopifnot(length(grep("'u0' must be a data.frame",
                      res[[1]]$message)) > 0)

u0 <- structure(list(S  = c(0, 1, 2, 3, 4, 5),
                     I  = c(0, 0, 0, 0, 0, 0)),
                .Names = c("S", "I"),
                row.names = c(NA, -6L), class = "data.frame")

## Check missing columns in u0
res <- tools::assertError(SISe(u0 = u0[, "I", drop = FALSE]))
stopifnot(length(grep("Missing columns in u0",
                      res[[1]]$message)) > 0)
res <- tools::assertError(SISe(u0 = u0[, "S", drop = FALSE]))
stopifnot(length(grep("Missing columns in u0",
                      res[[1]]$message)) > 0)

## Check phi
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep("a", nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("Invalid 'phi': must be numeric vector",
                      res[[1]]$message)) > 0)

res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = matrix(rep(1, nrow(u0))),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("Invalid 'phi': must be numeric vector",
                      res[[1]]$message)) > 0)

res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0) - 1),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(
    grep("Invalid 'phi': must be numeric vector with length 'nrow[(]u0[)]'",
                      res[[1]]$message)) > 0)

res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(-1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(
    grep("Invalid 'phi': must be numeric vector with non-negative values",
         res[[1]]$message)) > 0)

## Check missing upsilon
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'upsilon' is missing",
                      res[[1]]$message)) > 0)

## Check missing gamma
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'gamma' is missing",
                      res[[1]]$message)) > 0)

## Check missing alpha
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'alpha' is missing",
                      res[[1]]$message)) > 0)

## Check missing beta_t1
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'beta_t1' is missing",
                      res[[1]]$message)) > 0)

## Check missing beta_t2
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'beta_t2' is missing",
                      res[[1]]$message)) > 0)

## Check missing beta_t3
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'beta_t3' is missing",
                      res[[1]]$message)) > 0)

## Check missing beta_t4
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'beta_t4' is missing",
                      res[[1]]$message)) > 0)

## Check missing end_t1
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'end_t1' is missing",
                      res[[1]]$message)) > 0)

## Check missing end_t2
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'end_t2' is missing",
                      res[[1]]$message)) > 0)

## Check missing end_t3
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'end_t3' is missing",
                      res[[1]]$message)) > 0)

## Check missing end_t4
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               epsilon = 0.000011))
stopifnot(length(grep("'end_t4' is missing",
                      res[[1]]$message)) > 0)

## Check missing epsilon
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365))
stopifnot(length(grep("'epsilon' is missing",
                      res[[1]]$message)) > 0)

## Check non-numeric upsilon
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = "0.0357",
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'upsilon' must be numeric",
                      res[[1]]$message)) > 0)

## Check non-numeric gamma
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = "0.1",
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'gamma' must be numeric",
                      res[[1]]$message)) > 0)

## Check non-numeric alpha
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = "1.0",
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'alpha' must be numeric",
                      res[[1]]$message)) > 0)

## Check non-numeric beta_t1
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = "0.19",
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'beta_t1' must be numeric",
                      res[[1]]$message)) > 0)

## Check non-numeric beta_t2
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = "0.085",
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'beta_t2' must be numeric",
                      res[[1]]$message)) > 0)

## Check non-numeric beta_t3
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = "0.075",
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'beta_t3' must be numeric",
                      res[[1]]$message)) > 0)

## Check non-numeric beta_t4
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = "0.185",
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'beta_t4' must be numeric",
                      res[[1]]$message)) > 0)

## Check non-integer end_t1
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = "91",
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'end_t1' must be integer",
                      res[[1]]$message)) > 0)

res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91.5,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'end_t1' must be integer",
                      res[[1]]$message)) > 0)

## Check non-integer end_t2
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = "182",
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'end_t2' must be integer",
                      res[[1]]$message)) > 0)

res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182.5,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'end_t2' must be integer",
                      res[[1]]$message)) > 0)

## Check non-integer end_t3
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = "273",
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'end_t3' must be integer",
                      res[[1]]$message)) > 0)

res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273.5,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'end_t3' must be integer",
                      res[[1]]$message)) > 0)

## Check non-integer end_t4
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = "365",
                               epsilon = 0.000011))
stopifnot(length(grep("'end_t4' must be integer",
                      res[[1]]$message)) > 0)

res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365.5,
                               epsilon = 0.000011))
stopifnot(length(grep("'end_t4' must be integer",
                      res[[1]]$message)) > 0)

## Check non-numeric epsilon
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = "0.000011"))
stopifnot(length(grep("'epsilon' must be numeric",
                      res[[1]]$message)) > 0)

## Check that length of upsilon equals 1
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = c(0.0357, 0.0357),
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'upsilon' must be of length 1",
                      res[[1]]$message)) > 0)

## Check that length of gamma equals 1
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = c(0.1, 0.1),
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'gamma' must be of length 1",
                      res[[1]]$message)) > 0)

## Check that length of alpha equals 1
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = c(1.0, 1.0),
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'alpha' must be of length 1",
                      res[[1]]$message)) > 0)

## Check that length of beta_t1 equals 1
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = c(0.19, 0.19),
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'beta_t1' must be of length 1",
                      res[[1]]$message)) > 0)

## Check that length of beta_t2 equals 1
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = c(0.085, 0.085),
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'beta_t2' must be of length 1",
                      res[[1]]$message)) > 0)

## Check that length of beta_t3 equals 1
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = c(0.075, 0.075),
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'beta_t3' must be of length 1",
                      res[[1]]$message)) > 0)

## Check that length of beta_t4 equals 1
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = c(0.185, 0.185),
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'beta_t4' must be of length 1",
                      res[[1]]$message)) > 0)

## Check that length of end_t1 equals 1 or nrow(u0)
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = c(91, 91),
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'end_t1' must be of length 1 or 'nrow\\(u0\\)'",
                      res[[1]]$message)) > 0)

## Check that length of end_t2 equals 1 or nrow(u0)
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = c(182, 182),
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'end_t2' must be of length 1 or 'nrow\\(u0\\)'",
                      res[[1]]$message)) > 0)

## Check that length of end_t3 equals 1 or nrow(u0)
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = c(273, 273),
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'end_t3' must be of length 1 or 'nrow\\(u0\\)'",
                      res[[1]]$message)) > 0)

## Check that length of end_t4 equals 1 or nrow(u0)
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = c(365, 365),
                               epsilon = 0.000011))
stopifnot(length(grep("'end_t4' must be of length 1 or 'nrow\\(u0\\)'",
                      res[[1]]$message)) > 0)

## Check that length of epsilon equals 1
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = c(0.000011, 0.000011)))
stopifnot(length(grep("'epsilon' must be of length 1",
                      res[[1]]$message)) > 0)

## Check interval endpoints
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = -1,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'end_t1' must be greater than or equal to '0'",
                      res[[1]]$message)) > 0)

res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 18,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'end_t1' must be less than 'end_t2'",
                      res[[1]]$message)) > 0)

res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 173,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'end_t2' must be less than 'end_t3'",
                      res[[1]]$message)) > 0)

res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 365,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'end_t3' must be less than '364'",
                      res[[1]]$message)) > 0)

res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = -1,
                               epsilon = 0.000011))
stopifnot(length(grep("'end_t4' must be greater than or equal to '0'",
                      res[[1]]$message)) > 0)

res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 366,
                               epsilon = 0.000011))
stopifnot(length(grep("'end_t4' must be less than or equal to '365'",
                      res[[1]]$message)) > 0)

res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 4:9,
                               end_t2  = 5:10,
                               end_t3  = c(8:12, 16),
                               end_t4  = c(2, 11:15),
                               epsilon = 0.000011))
stopifnot(length(grep(
    "'end_t4' must be less than 'end_t1' or greater than 'end_t3'",
    res[[1]]$message)) > 0)

## Check 'suscpetible' and 'infected' methods
model <- SISe(u0      = u0,
              tspan   = seq_len(10) - 1,
              events  = NULL,
              phi     = rep(0, nrow(u0)),
              upsilon = 0.0357,
              gamma   = 0.1,
              alpha   = 1.0,
              beta_t1 = 0.19,
              beta_t2 = 0.085,
              beta_t3 = 0.075,
              beta_t4 = 0.185,
              end_t1  = 91,
              end_t2  = 182,
              end_t3  = 273,
              end_t4  = 365,
              epsilon = 0.000011)

result <- run(model)

S_expected <- structure(c(0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L, 0L,
                          1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L,
                          2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L,
                          3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L,
                          4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L),
                        .Dim = c(6L, 10L))

S_observed <- susceptible(result)

stopifnot(identical(S_observed, S_expected))

I_expected <- structure(c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
                        .Dim = c(6L, 10L))

I_observed <- infected(result)

stopifnot(identical(I_observed, I_expected))

## Check SISe plot method
pdf_file <- tempfile(fileext = ".pdf")
pdf(pdf_file)
plot(result)
dev.off()
stopifnot(file.exists(pdf_file))
unlink(pdf_file)

## Check that C SISe run function fails for misspecified SISe model
res <- .Call("SISe_run", NULL, NULL, NULL, PACKAGE = "SimInf")
stopifnot(identical(res$error, -10L))

res <- .Call("SISe_run", "SISe", NULL, NULL, PACKAGE = "SimInf")
stopifnot(identical(res$error, -10L))

## Check error non-finite v
model <- SISe(u0      = u0,
              tspan   = seq_len(10) - 1,
              events  = NULL,
              phi     = rep(1, nrow(u0)),
              upsilon = 0.0357,
              gamma   = 0.1,
              alpha   = 1.0,
              beta_t1 = 0.19,
              beta_t2 = 0.085,
              beta_t3 = 0.075,
              beta_t4 = 0.185,
              end_t1  = 91,
              end_t2  = 182,
              end_t3  = 273,
              end_t4  = 365,
              epsilon = 0.000011)
model@gdata <- rep(Inf, length(model@gdata))
res <- tools::assertError(run(model))
stopifnot(length(grep("The continuous state 'v' is not finite.",
                      res[[1]]$message)) > 0)
