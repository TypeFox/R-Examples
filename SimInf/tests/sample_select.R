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

## Test that sample_select in events.c works

## 2 Nodes
## 3 Age categories
## 2 Disease-states: Susceptible & Infected
##
## One individual start in susceptible state in node = 2, with a zero
## probability of becoming infected.
##
## At t = 1, two individuals are moved to node = 1. This should fail.
u0 <- structure(list(S_1 = c(0, 1),
                     I_1 = c(0, 0),
                     S_2 = c(0, 0),
                     I_2 = c(0, 0),
                     S_3 = c(0, 0),
                     I_3 = c(0, 0)),
                .Names = c("S_1", "I_1", "S_2", "I_2", "S_3", "I_3"),
                row.names = c(NA, -2L),
                class = "data.frame")

events <- structure(list(event      = 3,
                         time       = 1,
                         node       = 2,
                         dest       = 1,
                         n          = 2,
                         proportion = 1,
                         select     = 4,
                         shift      = 0),
                    .Names = c("event", "time", "node", "dest",
                        "n", "proportion", "select", "shift"),
                    row.names = c(NA, -1L), class = "data.frame")

model <- SISe3(u0        = u0,
               tspan     = 0:10,
               events    = events,
               phi       = rep(0, 2),
               upsilon_1 = 0,
               upsilon_2 = 0,
               upsilon_3 = 0,
               gamma_1   = 1,
               gamma_2   = 1,
               gamma_3   = 1,
               alpha     = 0,
               beta_t1   = 1,
               beta_t2   = 1,
               beta_t3   = 1,
               beta_t4   = 1,
               end_t1    = 91,
               end_t2    = 182,
               end_t3    = 273,
               end_t4    = 365,
               epsilon   = 0)

res <- tools::assertError(run(model, threads = 1))
stopifnot(length(grep("Unable to sample individuals for event.",
                      res[[1]]$message)) > 0)

if (SimInf:::have_openmp()) {
    res <- tools::assertError(run(model, threads = 2))
    stopifnot(length(grep("Unable to sample individuals for event.",
                          res[[1]]$message)) > 0)
}

## 2 Nodes
## 3 Age categories
## 2 Disease-states: Susceptible & Infected
##
## One individual start in susceptible state in node = 2, with a zero
## probability of becoming infected.
##
## At t = 1, -1 individuals are moved to node = 1. This should fail.
u0 <- structure(list(S_1 = c(0, 1),
                     I_1 = c(0, 0),
                     S_2 = c(0, 0),
                     I_2 = c(0, 0),
                     S_3 = c(0, 0),
                     I_3 = c(0, 0)),
                .Names = c("S_1", "I_1", "S_2", "I_2", "S_3", "I_3"),
                row.names = c(NA, -2L), class = "data.frame")

events <- structure(list(event      = 3,
                         time       = 1,
                         node       = 2,
                         dest       = 1,
                         n          = -1,
                         proportion = 1,
                         select     = 4,
                         shift      = 0),
                    .Names = c("event", "time", "node", "dest",
                        "n", "proportion", "select", "shift"),
                    row.names = c(NA, -1L), class = "data.frame")

model <- SISe3(u0        = u0,
               tspan     = 0:10,
               events    = events,
               phi       = rep(0, 2),
               upsilon_1 = 0,
               upsilon_2 = 0,
               upsilon_3 = 0,
               gamma_1   = 1,
               gamma_2   = 1,
               gamma_3   = 1,
               alpha     = 0,
               beta_t1   = 1,
               beta_t2   = 1,
               beta_t3   = 1,
               beta_t4   = 1,
               end_t1    = 91,
               end_t2    = 182,
               end_t3    = 273,
               end_t4    = 365,
               epsilon   = 0)

res <- tools::assertError(run(model, threads = 1))
stopifnot(length(grep("Unable to sample individuals for event.",
                      res[[1]]$message)) > 0)

if (SimInf:::have_openmp()) {
    res <- tools::assertError(run(model, threads = 2))
    stopifnot(length(grep("Unable to sample individuals for event.",
                          res[[1]]$message)) > 0)
}

## 2 Nodes
## 3 Age categories
## 2 Disease-states: Susceptible & Infected
##
## One individual start in susceptible state in node = 2, with a zero
## probability of becoming infected.
##
## At t = 1, a proportion of 10 individuals are moved to node = 1.
## This should fail.
u0 <- structure(list(S_1 = c(0, 1),
                     I_1 = c(0, 0),
                     S_2 = c(0, 0),
                     I_2 = c(0, 0),
                     S_3 = c(0, 0),
                     I_3 = c(0, 0)),
                .Names = c("S_1", "I_1", "S_2", "I_2", "S_3", "I_3"),
                row.names = c(NA, -2L), class = "data.frame")

events <- structure(list(event      = 3,
                         time       = 1,
                         node       = 2,
                         dest       = 1,
                         n          = 0,
                         proportion = 10,
                         select     = 4,
                         shift      = 0),
                    .Names = c("event", "time", "node", "dest",
                        "n", "proportion", "select", "shift"),
                    row.names = c(NA, -1L), class = "data.frame")

## We should not be able to create model with prop = 10
res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = 0:10,
                                events    = events,
                                phi       = rep(0, 2),
                                upsilon_1 = 0,
                                upsilon_2 = 0,
                                upsilon_3 = 0,
                                gamma_1   = 1,
                                gamma_2   = 1,
                                gamma_3   = 1,
                                alpha     = 0,
                                beta_t1   = 1,
                                beta_t2   = 1,
                                beta_t3   = 1,
                                beta_t4   = 1,
                                end_t1    = 91,
                                end_t2    = 182,
                                end_t3    = 273,
                                end_t4    = 365,
                                epsilon   = 0))
stopifnot(length(grep("prop must be in the range 0 <= prop <= 1",
                      res[[1]]$message)) > 0)

## Replace proportion = 10 to proportion = 1
events$proportion <- 1

model <- SISe3(u0        = u0,
               tspan     = 0:10,
               events    = events,
               phi       = rep(0, 2),
               upsilon_1 = 0,
               upsilon_2 = 0,
               upsilon_3 = 0,
               gamma_1   = 1,
               gamma_2   = 1,
               gamma_3   = 1,
               alpha     = 0,
               beta_t1   = 1,
               beta_t2   = 1,
               beta_t3   = 1,
               beta_t4   = 1,
               end_t1    = 91,
               end_t2    = 182,
               end_t3    = 273,
               end_t4    = 365,
               epsilon   = 0)

## Replace proportion = 1 with proportion = 10
model@events@proportion <- 10

res <- .Call("SISe3_run", model, 1L, NULL, PACKAGE = "SimInf")
res <- tools::assertError(SimInf:::siminf_error(res$error))
stopifnot(length(grep("Unable to sample individuals for event.",
                      res[[1]]$message)) > 0)

if (SimInf:::have_openmp()) {
    res <- .Call("SISe3_run", model, 2L, NULL, PACKAGE = "SimInf")
    res <- tools::assertError(SimInf:::siminf_error(res$error))
    stopifnot(length(grep("Unable to sample individuals for event.",
                          res[[1]]$message)) > 0)
}

## 2 Nodes
## 3 Age categories
## 2 Disease-states: Susceptible & Infected
##
## One individual start in susceptible state in node = 2, with a zero
## probability of becoming infected.
##
## At t = 1, a proportion of -1 individuals are moved to node = 1.
## This should fail.
u0 <- structure(list(S_1 = c(0, 1),
                     I_1 = c(0, 0),
                     S_2 = c(0, 0),
                     I_2 = c(0, 0),
                     S_3 = c(0, 0),
                     I_3 = c(0, 0)),
                .Names = c("S_1", "I_1", "S_2", "I_2", "S_3", "I_3"),
                row.names = c(NA, -2L), class = "data.frame")

events <- structure(list(event      = 3,
                         time       = 1,
                         node       = 2,
                         dest       = 1,
                         n          = 0,
                         proportion = -1,
                         select     = 4,
                         shift      = 0),
                    .Names = c("event", "time", "node", "dest",
                        "n", "proportion", "select", "shift"),
                    row.names = c(NA, -1L), class = "data.frame")

## We should not be able to create model with proportion = -1
res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = 0:10,
                                events    = events,
                                phi       = rep(0, 2),
                                upsilon_1 = 0,
                                upsilon_2 = 0,
                                upsilon_3 = 0,
                                gamma_1   = 1,
                                gamma_2   = 1,
                                gamma_3   = 1,
                                alpha     = 0,
                                beta_t1   = 1,
                                beta_t2   = 1,
                                beta_t3   = 1,
                                beta_t4   = 1,
                                end_t1    = 91,
                                end_t2    = 182,
                                end_t3    = 273,
                                end_t4    = 365,
                                epsilon   = 0))
stopifnot(length(grep("prop must be in the range 0 <= prop <= 1",
                      res[[1]]$message)) > 0)

## Replace proportion = -1 with proportion = 0
events$proportion <- 0

model <- SISe3(u0        = u0,
               tspan     = 0:10,
               events    = events,
               phi       = rep(0, 2),
               upsilon_1 = 0,
               upsilon_2 = 0,
               upsilon_3 = 0,
               gamma_1   = 1,
               gamma_2   = 1,
               gamma_3   = 1,
               alpha     = 0,
               beta_t1   = 1,
               beta_t2   = 1,
               beta_t3   = 1,
               beta_t4   = 1,
               end_t1    = 91,
               end_t2    = 182,
               end_t3    = 273,
               end_t4    = 365,
               epsilon   = 0)

## Replace proportion = 0 with proportion = -1
model@events@proportion <- -1

res <- .Call("SISe3_run", model, 1L, NULL, PACKAGE = "SimInf")
res <- tools::assertError(SimInf:::siminf_error(res$error))
stopifnot(length(grep("Unable to sample individuals for event.",
                      res[[1]]$message)) > 0)

if (SimInf:::have_openmp()) {
    res <- .Call("SISe3_run", model, 2L, NULL, PACKAGE = "SimInf")
    res <- tools::assertError(SimInf:::siminf_error(res$error))
    stopifnot(length(grep("Unable to sample individuals for event.",
                          res[[1]]$message)) > 0)
}

## 2 Nodes
## 3 Age categories
## 2 Disease-states: Susceptible & Infected
##
## One individual start in susceptible state in node = 2, with a zero
## probability of becoming infected.
##
## At t = 1, a proportion of 0 individuals are moved to node = 1.
u0 <- structure(list(S_1 = c(0, 1),
                     I_1 = c(0, 0),
                     S_2 = c(0, 0),
                     I_2 = c(0, 0),
                     S_3 = c(0, 0),
                     I_3  = c(0, 0)),
                .Names = c("S_1", "I_1", "S_2", "I_2", "S_3", "I_3"),
                row.names = c(NA, -2L), class = "data.frame")

events <- structure(list(event      = 3,
                         time       = 1,
                         node       = 2,
                         dest       = 1,
                         n          = 0,
                         proportion = 0,
                         select     = 4,
                         shift      = 0),
                    .Names = c("event", "time", "node", "dest",
                        "n", "proportion", "select", "shift"),
                    row.names = c(NA, -1L), class = "data.frame")

model <- SISe3(u0        = u0,
               tspan     = 0:2,
               events    = events,
               phi       = rep(0, 2),
               upsilon_1 = 0,
               upsilon_2 = 0,
               upsilon_3 = 0,
               gamma_1   = 1,
               gamma_2   = 1,
               gamma_3   = 1,
               alpha     = 0,
               beta_t1   = 1,
               beta_t2   = 1,
               beta_t3   = 1,
               beta_t4   = 1,
               end_t1    = 91,
               end_t2    = 182,
               end_t3    = 273,
               end_t4    = 365,
               epsilon   = 0)

U <- structure(c(0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 1L, 0L, 0L, 0L, 0L, 0L), .Dim = c(12L, 3L))

result <- run(model, threads = 1)
stopifnot(identical(model@G, result@G))
stopifnot(identical(model@S, result@S))
stopifnot(identical(result@U, U))
stopifnot(identical(model@ldata, result@ldata))
stopifnot(identical(model@sd, result@sd))
stopifnot(identical(model@tspan, result@tspan))
stopifnot(identical(model@u0, result@u0))
stopifnot(identical(model@events, result@events))

if (SimInf:::have_openmp()) {
    result_omp <- run(model, threads = 2)
    stopifnot(identical(model@G, result_omp@G))
    stopifnot(identical(model@S, result_omp@S))
    stopifnot(identical(result_omp@U, U))
    stopifnot(identical(model@ldata, result_omp@ldata))
    stopifnot(identical(model@sd, result_omp@sd))
    stopifnot(identical(model@tspan, result_omp@tspan))
    stopifnot(identical(model@u0, result_omp@u0))
    stopifnot(identical(model@events, result_omp@events))
}

## 2 Nodes
## 3 Age categories
## 2 Disease-states: Susceptible & Infected
##
## One individual start in susceptible state in node = 2, with a zero
## probability of becoming infected.
##
## At t = 1, proportion of all (1) individuals are moved to node = 1.
u0 <- structure(list(S_1 = c(0, 1),
                     I_1 = c(0, 0),
                     S_2 = c(0, 0),
                     I_2 = c(0, 0),
                     S_3  = c(0, 0),
                     I_3  = c(0, 0)),
                .Names = c("S_1", "I_1", "S_2", "I_2", "S_3", "I_3"),
                row.names = c(NA, -2L), class = "data.frame")

events <- structure(list(event      = 3,
                         time       = 1,
                         node       = 2,
                         dest       = 1,
                         n          = 1,
                         proportion = 0,
                         select     = 4,
                         shift      = 0),
                    .Names = c("event", "time", "node", "dest",
                        "n", "proportion", "select", "shift"),
                    row.names = c(NA, -1L), class = "data.frame")

model <- SISe3(u0        = u0,
               tspan     = 0:2,
               events    = events,
               phi       = rep(0, 2),
               upsilon_1 = 0,
               upsilon_2 = 0,
               upsilon_3 = 0,
               gamma_1   = 1,
               gamma_2   = 1,
               gamma_3   = 1,
               alpha     = 0,
               beta_t1   = 1,
               beta_t2   = 1,
               beta_t3   = 1,
               beta_t4   = 1,
               end_t1    = 91,
               end_t2    = 182,
               end_t3    = 273,
               end_t4    = 365,
               epsilon   = 0)

U <- structure(c(0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 1L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L), .Dim = c(12L, 3L))

result <- run(model, threads = 1)
stopifnot(identical(model@G, result@G))
stopifnot(identical(model@S, result@S))
stopifnot(identical(result@U, U))
stopifnot(identical(model@ldata, result@ldata))
stopifnot(identical(model@sd, result@sd))
stopifnot(identical(model@tspan, result@tspan))
stopifnot(identical(model@u0, result@u0))
stopifnot(identical(model@events, result@events))

if (SimInf:::have_openmp()) {
    result_omp <- run(model, threads = 2)
    stopifnot(identical(model@G, result_omp@G))
    stopifnot(identical(model@S, result_omp@S))
    stopifnot(identical(result_omp@U, U))
    stopifnot(identical(model@ldata, result_omp@ldata))
    stopifnot(identical(model@sd, result_omp@sd))
    stopifnot(identical(model@tspan, result_omp@tspan))
    stopifnot(identical(model@u0, result_omp@u0))
    stopifnot(identical(model@events, result_omp@events))
}

## 2 Nodes
## 3 Age categories
## 2 Disease-states: Susceptible & Infected
##
## Two individuals start in susceptible state in node = 2, with a zero
## probability of becoming infected.
##
## At t = 1, Nkind = 1, one individual is moved to node = 1.
u0 <- structure(list(S_1 = c(0, 2),
                     I_1 = c(0, 0),
                     S_2 = c(0, 0),
                     I_2 = c(0, 0),
                     S_3 = c(0, 0),
                     I_3  = c(0, 0)),
                .Names = c("S_1", "I_1", "S_2", "I_2", "S_3", "I_3"),
                row.names = c(NA, -2L), class = "data.frame")

events <- structure(list(event      = 3,
                         time       = 1,
                         node       = 2,
                         dest       = 1,
                         n          = 1,
                         proportion = 0,
                         select     = 4,
                         shift      = 0),
                    .Names = c("event", "time", "node", "dest",
                        "n", "proportion", "select", "shift"),
                    row.names = c(NA, -1L), class = "data.frame")

model <- SISe3(u0        = u0,
               tspan     = 0:2,
               events    = events,
               phi       = rep(0, 2),
               upsilon_1 = 0,
               upsilon_2 = 0,
               upsilon_3 = 0,
               gamma_1   = 1,
               gamma_2   = 1,
               gamma_3   = 1,
               alpha     = 0,
               beta_t1   = 1,
               beta_t2   = 1,
               beta_t3   = 1,
               beta_t4   = 1,
               end_t1    = 91,
               end_t2    = 182,
               end_t3    = 273,
               end_t4    = 365,
               epsilon   = 0)

U <- structure(c(0L, 0L, 0L, 0L, 0L, 0L, 2L, 0L, 0L, 0L, 0L, 0L, 1L,
                 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L,
                 0L, 1L, 0L, 0L, 0L, 0L, 0L), .Dim = c(12L, 3L))

result <- run(model, threads = 1)
stopifnot(identical(model@G, result@G))
stopifnot(identical(model@S, result@S))
stopifnot(identical(result@U, U))
stopifnot(identical(model@ldata, result@ldata))
stopifnot(identical(model@sd, result@sd))
stopifnot(identical(model@tspan, result@tspan))
stopifnot(identical(model@u0, result@u0))
stopifnot(identical(model@events, result@events))

if (SimInf:::have_openmp()) {
    result_omp <- run(model, threads = 2)
    stopifnot(identical(model@G, result_omp@G))
    stopifnot(identical(model@S, result_omp@S))
    stopifnot(identical(result_omp@U, U))
    stopifnot(identical(model@ldata, result_omp@ldata))
    stopifnot(identical(model@sd, result_omp@sd))
    stopifnot(identical(model@tspan, result_omp@tspan))
    stopifnot(identical(model@u0, result_omp@u0))
    stopifnot(identical(model@events, result_omp@events))
}

## 2 Nodes
## 3 Age categories
## 2 Disease-states: Susceptible & Infected
##
## Two individuals start in susceptible state and 8 individuals in
## infected state in node = 2, with a zero probability of becoming
## infected.
##
## At t = 1, Nkind = 2, one individual is moved to node = 1.
u0 <- structure(list(S_1 = c(0, 2),
                     I_1 = c(0, 8),
                     S_2 = c(0, 0),
                     I_2 = c(0, 0),
                     S_3 = c(0, 0),
                     I_3 = c(0, 0)),
                .Names = c("S_1", "I_1", "S_2", "I_2", "S_3", "I_3"),
                row.names = c(NA, -2L), class = "data.frame")

events <- structure(list(event      = 3,
                         time       = 1,
                         node       = 2,
                         dest       = 1,
                         n          = 1,
                         proportion = 0,
                         select     = 4,
                         shift      = 0),
                    .Names = c("event", "time", "node", "dest",
                        "n", "proportion", "select", "shift"),
                    row.names = c(NA, -1L), class = "data.frame")

model <- SISe3(u0        = u0,
               tspan     = 0:2,
               events    = events,
               phi       = rep(0, 2),
               upsilon_1 = 0,
               upsilon_2 = 0,
               upsilon_3 = 0,
               gamma_1   = 0,
               gamma_2   = 0,
               gamma_3   = 0,
               alpha     = 0,
               beta_t1   = 1,
               beta_t2   = 1,
               beta_t3   = 1,
               beta_t4   = 1,
               end_t1    = 91,
               end_t2    = 182,
               end_t3    = 273,
               end_t4    = 365,
               epsilon   = 0)

U <- structure(c(0L, 0L, 0L, 0L, 0L, 0L, 2L, 8L, 0L, 0L, 0L, 0L, 0L,
                 1L, 0L, 0L, 0L, 0L, 2L, 7L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L,
                 0L, 2L, 7L, 0L, 0L, 0L, 0L), .Dim = c(12L, 3L))

result <- run(model, threads = 1, seed = 123L)
stopifnot(identical(model@G, result@G))
stopifnot(identical(model@S, result@S))
stopifnot(identical(result@U, U))
stopifnot(identical(model@ldata, result@ldata))
stopifnot(identical(model@sd, result@sd))
stopifnot(identical(model@tspan, result@tspan))
stopifnot(identical(model@u0, result@u0))
stopifnot(identical(model@events, result@events))

if (SimInf:::have_openmp()) {
    result_omp <- run(model, threads = 2, seed = 123L)
    stopifnot(identical(model@G, result_omp@G))
    stopifnot(identical(model@S, result_omp@S))
    stopifnot(identical(result_omp@U, U))
    stopifnot(identical(model@ldata, result_omp@ldata))
    stopifnot(identical(model@sd, result_omp@sd))
    stopifnot(identical(model@tspan, result_omp@tspan))
    stopifnot(identical(model@u0, result_omp@u0))
    stopifnot(identical(model@events, result_omp@events))
}

## 2 Nodes
## 3 Age categories
## 2 Disease-states: Susceptible & Infected
##
## Two individuals start in susceptible state in node = 2, with a zero
## probability of becoming infected.
##
## At t = 1, one individual is moved to node = 1 using a select matrix
## that can only select a susceptible from S_1.
u0 <- structure(list(S_1 = c(0, 2),
                     I_1 = c(0, 0),
                     S_2 = c(0, 0),
                     I_2 = c(0, 0),
                     S_3 = c(0, 0),
                     I_3 = c(0, 0)),
                .Names = c("S_1", "I_1", "S_2", "I_2", "S_3", "I_3"),
                row.names = c(NA, -2L),
                class = "data.frame")

events <- structure(list(event      = 3,
                         time       = 1,
                         node       = 2,
                         dest       = 1,
                         n          = 1,
                         proportion = 0,
                         select     = 1,
                         shift      = 0),
                    .Names = c("event", "time", "node", "dest",
                        "n", "proportion", "select", "shift"),
                    row.names = c(NA, -1L), class = "data.frame")

model <- SISe3(u0        = u0,
               tspan     = 0:10,
               events    = events,
               phi       = rep(0, 2),
               upsilon_1 = 0,
               upsilon_2 = 0,
               upsilon_3 = 0,
               gamma_1   = 1,
               gamma_2   = 1,
               gamma_3   = 1,
               alpha     = 0,
               beta_t1   = 1,
               beta_t2   = 1,
               beta_t3   = 1,
               beta_t4   = 1,
               end_t1    = 91,
               end_t2    = 182,
               end_t3    = 273,
               end_t4    = 365,
               epsilon   = 0)

U <- structure(
    c(0L, 0L, 0L, 0L, 0L, 0L, 2L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L,
      0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 0L,
      0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L,
      1L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L,
      0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 0L,
      0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L,
      1L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L,
      0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 0L,
      0L, 0L, 0L, 0L),
    .Dim = c(12L, 11L))

res <- run(model, threads = 1)
stopifnot(identical(res@U, U))

if (SimInf:::have_openmp()) {
    res <- run(model, threads = 2)
    stopifnot(identical(res@U, U))
}

## 2 Nodes
## 3 Age categories
## 2 Disease-states: Susceptible & Infected
##
## 6 individuals start in the susceptible compartments in node = 2,
## with a zero probability of becoming infected.
##
## At t = 1, three individuals are moved to node = 1 using a select
## matrix that can select from the susceptible S_1, S_2 and S_3.
u0 <- structure(list(S_1 = c(0, 2),
                     I_1 = c(0, 0),
                     S_2 = c(0, 2),
                     I_2 = c(0, 0),
                     S_3 = c(0, 2),
                     I_3 = c(0, 0)),
                .Names = c("S_1", "I_1", "S_2", "I_2", "S_3", "I_3"),
                row.names = c(NA, -2L),
                class = "data.frame")

events <- structure(list(event      = 3,
                         time       = 1,
                         node       = 2,
                         dest       = 1,
                         n          = 3,
                         proportion = 0,
                         select     = 1,
                         shift      = 0),
                    .Names = c("event", "time", "node", "dest",
                        "n", "proportion", "select", "shift"),
                    row.names = c(NA, -1L), class = "data.frame")

model <- SISe3(u0        = u0,
               tspan     = 0:10,
               events    = events,
               phi       = rep(0, 2),
               upsilon_1 = 0,
               upsilon_2 = 0,
               upsilon_3 = 0,
               gamma_1   = 1,
               gamma_2   = 1,
               gamma_3   = 1,
               alpha     = 0,
               beta_t1   = 1,
               beta_t2   = 1,
               beta_t3   = 1,
               beta_t4   = 1,
               end_t1    = 91,
               end_t2    = 182,
               end_t3    = 273,
               end_t4    = 365,
               epsilon   = 0)

## Add S_2 and S_3 to first column in select matrix
model@events@E[3, 1] <- 1
model@events@E[5, 1] <- 1

S_expected <- structure(c(0L, 6L, 3L, 3L, 3L, 3L, 3L, 3L, 3L,
                          3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L,
                          3L, 3L, 3L, 3L),
                        .Dim = c(2L, 11L))

res <- run(model, threads = 1)
stopifnot(identical(susceptible(res), S_expected))

if (SimInf:::have_openmp()) {
    res <- run(model, threads = 2)
    stopifnot(identical(susceptible(res), S_expected))
}
