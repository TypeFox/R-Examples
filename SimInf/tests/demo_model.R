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

res <- tools::assertError(demo_model(nodes = 1:2))
stopifnot(length(grep("Length of 'nodes' must be one.",
                      res[[1]]$message)) > 0)

res <- tools::assertError(demo_model(nodes = "1"))
stopifnot(length(grep("'nodes' must be numeric.",
                      res[[1]]$message)) > 0)

res <- tools::assertError(demo_model(nodes = 1.1))
stopifnot(length(grep("'nodes' must be integer",
                      res[[1]]$message)) > 0)

res <- tools::assertError(demo_model(nodes = 3, model = "SISe_sp"))
stopifnot(length(grep("'sqrt[(]nodes[)]' must be integer",
                      res[[1]]$message)) > 0)

res <- tools::assertError(demo_model(nodes = 0))
stopifnot(length(grep("'nodes' must be >= 1",
                      res[[1]]$message)) > 0)

res <- tools::assertError(demo_model(nodes = -1))
stopifnot(length(grep("'nodes' must be >= 1",
                      res[[1]]$message)) > 0)

res <- tools::assertError(demo_model(days = 1:2))
stopifnot(length(grep("Length of 'days' must be one.",
                      res[[1]]$message)) > 0)

res <- tools::assertError(demo_model(days = "1"))
stopifnot(length(grep("'days' must be numeric.",
                      res[[1]]$message)) > 0)

res <- tools::assertError(demo_model(days = 1.1))
stopifnot(length(grep("'days' must be integer",
                      res[[1]]$message)) > 0)

res <- tools::assertError(demo_model(days = 0))
stopifnot(length(grep("'days' must be >= 1",
                      res[[1]]$message)) > 0)

res <- tools::assertError(demo_model(days = -1))
stopifnot(length(grep("'days' must be >= 1",
                      res[[1]]$message)) > 0)

res <- tools::assertError(demo_model(model = "demo"))
stopifnot(length(grep("'arg' should be one of",
                      res[[1]]$message)) > 0)

stopifnot(is(demo_model(model = "SISe"), "SISe"))
stopifnot(is(demo_model(model = "SISe3"), "SISe3"))
stopifnot(is(demo_model(model = "SISe_sp"), "SISe_sp"))
