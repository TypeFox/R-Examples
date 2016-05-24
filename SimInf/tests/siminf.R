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

## 1 Node
## 1 Age category
## 2 Disease-states: Susceptible & Infected
##
## The individual start in the susceptible state, with a probability
## of becoming infected.
u0 <- structure(list(S = 1, I = 0),
                .Names = c("S", "I"),
                row.names = c(NA, -1L),
                class = "data.frame")

model <- SISe(u0      = u0,
              tspan   = 0:1000,
              events  = NULL,
              phi     = 1,
              upsilon = 1,
              gamma   = 1,
              alpha   = 1,
              beta_t1 = 1,
              beta_t2 = 1,
              beta_t3 = 1,
              beta_t4 = 1,
              end_t1  = 91,
              end_t2  = 182,
              end_t3  = 273,
              end_t4  = 365,
              epsilon = 1)

result <- run(model, threads = 1, seed = 123L)
stopifnot(identical(model@G, result@G))
stopifnot(identical(model@S, result@S))
stopifnot(identical(sum(result@U), 1001L))
stopifnot(any(result@U[1,]))
stopifnot(any(result@U[2,]))
stopifnot(identical(model@ldata, result@ldata))
stopifnot(identical(model@sd, result@sd))
stopifnot(identical(model@tspan, result@tspan))
stopifnot(identical(model@u0, result@u0))
stopifnot(identical(model@events, result@events))

if (SimInf:::have_openmp()) {
    result_omp <- run(model, threads = 2, seed = 123L)
    stopifnot(identical(model@G, result_omp@G))
    stopifnot(identical(model@S, result_omp@S))
    stopifnot(identical(sum(result_omp@U), 1001L))
    stopifnot(any(result_omp@U[1,]))
    stopifnot(any(result_omp@U[2,]))
    stopifnot(identical(model@ldata, result_omp@ldata))
    stopifnot(identical(model@sd, result_omp@sd))
    stopifnot(identical(model@tspan, result_omp@tspan))
    stopifnot(identical(model@u0, result_omp@u0))
    stopifnot(identical(model@events, result_omp@events))
}

## 6 Nodes
## 3 Age categories
## 2 Disease-states: Susceptible & Infected
##
## All individuals start in susceptible state, with a probability of
## becoming infected.
##
## At t = 1, all individuals are moved to node = 1.
u0 <- structure(list(S_1 = c(0, 1, 2, 3, 4, 5),
                     I_1 = c(0, 0, 0, 0, 0, 0),
                     S_2 = c(0, 1, 2, 3, 4, 5),
                     I_2 = c(0, 0, 0, 0, 0, 0),
                     S_3 = c(0, 1, 2, 3, 4, 5),
                     I_3 = c(0, 0, 0, 0, 0, 0)),
                .Names = c("S_1", "I_1", "S_2", "I_2", "S_3", "I_3"),
                row.names = c(NA, -6L),
                class = "data.frame")

events <- structure(list(
    event      = c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3),
    time       = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    node       = c(2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6),
    dest       = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    n          = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
    proportion = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    select     = c(4, 5, 6, 4, 5, 6, 4, 5, 6, 4, 5, 6, 4, 5, 6),
    shift      = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)),
    .Names = c("event", "time", "node", "dest",
               "n", "proportion", "select", "shift"),
    row.names = c(NA, -15L), class = "data.frame")

model <- SISe3(u0        = u0,
               tspan     = 0:10,
               events    = events,
               phi       = rep(1, 6),
               upsilon_1 = 1,
               upsilon_2 = 1,
               upsilon_3 = 1,
               gamma_1   = 1,
               gamma_2   = 1,
               gamma_3   = 1,
               alpha     = 1,
               beta_t1   = 1,
               beta_t2   = 1,
               beta_t3   = 1,
               beta_t4   = 1,
               end_t1    = 91,
               end_t2    = 182,
               end_t3    = 273,
               end_t4    = 365,
               epsilon   = 1)

result <- run(model, threads = 1, seed = 123L)
stopifnot(identical(model@G, result@G))
stopifnot(identical(model@S, result@S))
stopifnot(all(apply(result@U[1:6,], 1, any)))
stopifnot(identical(sum(result@U[1:6,11]), 45L))
stopifnot(identical(sum(result@U[,1]), 45L))
stopifnot(identical(model@ldata, result@ldata))
stopifnot(identical(model@sd, result@sd))
stopifnot(identical(model@tspan, result@tspan))
stopifnot(identical(model@u0, result@u0))
stopifnot(identical(model@events, result@events))

if (SimInf:::have_openmp()) {
    result_omp <- run(model, threads = 2, seed = 123L)
    stopifnot(identical(model@G, result_omp@G))
    stopifnot(identical(model@S, result_omp@S))
    stopifnot(all(apply(result_omp@U[1:6,], 1, any)))
    stopifnot(identical(sum(result_omp@U[1:6,11]), 45L))
    stopifnot(identical(sum(result_omp@U[,1]), 45L))
    stopifnot(identical(model@ldata, result_omp@ldata))
    stopifnot(identical(model@sd, result_omp@sd))
    stopifnot(identical(model@tspan, result_omp@tspan))
    stopifnot(identical(model@u0, result_omp@u0))
    stopifnot(identical(model@events, result_omp@events))
}

## 6 Nodes
## 3 Age categories
## 2 Disease-states: Susceptible & Infected
##
## All individuals start in susceptible state, with a probability of
## becoming infected.
##
## At t = 1, all individuals in age category 1 and 2 are moved to node
## = 1 and age.
u0 <- structure(list(S_1 = c(0, 1, 2, 3, 4, 5),
                     I_1 = c(0, 0, 0, 0, 0, 0),
                     S_2 = c(0, 1, 2, 3, 4, 5),
                     I_2 = c(0, 0, 0, 0, 0, 0),
                     S_3 = c(0, 1, 2, 3, 4, 5),
                     I_3 = c(0, 0, 0, 0, 0, 0)),
                .Names = c("S_1", "I_1", "S_2", "I_2", "S_3", "I_3"),
                row.names = c(NA, -6L),
                class = "data.frame")

events <- structure(list(
    event      = c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3),
    time       = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    node       = c(2, 2, 3, 3, 4, 4, 5, 5, 6, 6),
    dest       = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    n          = c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5),
    proportion = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    select     = c(4, 5, 4, 5, 4, 5, 4, 5, 4, 5),
    shift      = c(1, 2, 1, 2, 1, 2, 1, 2, 1, 2)),
    .Names = c("event", "time", "node", "dest",
               "n", "proportion", "select", "shift"),
    row.names = c(NA, -15L), class = "data.frame")

model <- SISe3(u0        = u0,
               tspan     = 0:10,
               events    = events,
               phi       = rep(1, 6),
               upsilon_1 = 1,
               upsilon_2 = 1,
               upsilon_3 = 1,
               gamma_1   = 1,
               gamma_2   = 1,
               gamma_3   = 1,
               alpha     = 1,
               beta_t1   = 1,
               beta_t2   = 1,
               beta_t3   = 1,
               beta_t4   = 1,
               end_t1    = 91,
               end_t2    = 182,
               end_t3    = 273,
               end_t4    = 365,
               epsilon   = 1)

result <- run(model, threads = 1, seed = 123L)
stopifnot(identical(model@G, result@G))
stopifnot(identical(model@S, result@S))
stopifnot(identical(
    susceptible(result, age = 2, i = 1) + infected(result, age = 2, i = 1),
    structure(c(0L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L),
              .Dim = c(1L, 11L))))
stopifnot(identical(
    susceptible(result, age = 3, i = 1) + infected(result, age = 3, i = 1),
    structure(c(0L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L),
              .Dim = c(1L, 11L))))
stopifnot(identical(sum(result@U[,1]), 45L))
stopifnot(identical(model@ldata, result@ldata))
stopifnot(identical(model@sd, result@sd))
stopifnot(identical(model@tspan, result@tspan))
stopifnot(identical(model@u0, result@u0))
stopifnot(identical(model@events, result@events))

if (SimInf:::have_openmp()) {
    result_omp <- run(model, threads = 2, seed = 123L)
    stopifnot(identical(model@G, result_omp@G))
    stopifnot(identical(model@S, result_omp@S))
    stopifnot(identical(
        susceptible(result, age = 2, i = 1) + infected(result, age = 2, i = 1),
        structure(c(0L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L),
                  .Dim = c(1L, 11L))))
    stopifnot(identical(
        susceptible(result, age = 3, i = 1) + infected(result, age = 3, i = 1),
        structure(c(0L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L),
                  .Dim = c(1L, 11L))))
    stopifnot(identical(sum(result_omp@U[,1]), 45L))
    stopifnot(identical(model@ldata, result_omp@ldata))
    stopifnot(identical(model@sd, result_omp@sd))
    stopifnot(identical(model@tspan, result_omp@tspan))
    stopifnot(identical(model@u0, result_omp@u0))
    stopifnot(identical(model@events, result_omp@events))
}

## 6 Nodes
## 3 Age categories
## 2 Disease-states: Susceptible & Infected
##
## All individuals start in susceptible state, with a probability of
## becoming infected.
##
## No scheduled events
u0 <- structure(list(S_1 = c(0, 1, 2, 3, 4, 5),
                     I_1 = c(0, 0, 0, 0, 0, 0),
                     S_2 = c(0, 1, 2, 3, 4, 5),
                     I_2 = c(0, 0, 0, 0, 0, 0),
                     S_3 = c(0, 1, 2, 3, 4, 5),
                     I_3 = c(0, 0, 0, 0, 0, 0)),
                .Names = c("S_1", "I_1", "S_2", "I_2", "S_3", "I_3"),
                row.names = c(NA, -6L), class = "data.frame")

model <- SISe3(u0        = u0,
               tspan     = 0:10,
               events    = NULL,
               phi       = rep(1, 6),
               upsilon_1 = 1,
               upsilon_2 = 1,
               upsilon_3 = 1,
               gamma_1   = 1,
               gamma_2   = 1,
               gamma_3   = 1,
               alpha     = 1,
               beta_t1   = 1,
               beta_t2   = 1,
               beta_t3   = 1,
               beta_t4   = 1,
               end_t1    = 91,
               end_t2    = 182,
               end_t3    = 273,
               end_t4    = 365,
               epsilon   = 1)

result <- run(model, threads = 1, seed = 123L)
stopifnot(identical(model@G, result@G))
stopifnot(identical(model@S, result@S))
stopifnot(all(result@U[1:6,] == 0))
stopifnot(all(apply(result@U[seq(from=8, to=36, by=2),], 1, any)))
stopifnot(identical(sum(result@U[,11]), 45L))
stopifnot(identical(model@ldata, result@ldata))
stopifnot(identical(model@sd, result@sd))
stopifnot(identical(model@tspan, result@tspan))
stopifnot(identical(model@u0, result@u0))
stopifnot(identical(model@events, result@events))

if (SimInf:::have_openmp()) {
    result_omp <- run(model, threads = 2, seed = 123L)
    stopifnot(identical(model@G, result_omp@G))
    stopifnot(identical(model@S, result_omp@S))
    stopifnot(all(result_omp@U[1:6,] == 0))
    stopifnot(all(apply(result_omp@U[seq(from=8, to=36, by=2),], 1, any)))
    stopifnot(identical(sum(result_omp@U[,11]), 45L))
    stopifnot(identical(model@ldata, result_omp@ldata))
    stopifnot(identical(model@sd, result_omp@sd))
    stopifnot(identical(model@tspan, result_omp@tspan))
    stopifnot(identical(model@u0, result_omp@u0))
    stopifnot(identical(model@events, result_omp@events))
}

## 6 Nodes
## 3 Age categories
## 2 Disease-states: Susceptible & Infected
##
## All individuals start in susceptible state, with a zero probability
## of becoming infected.
##
## At t = 1, all individuals are moved to node = 1.
u0 <- structure(list(S_1 = c(0, 1, 2, 3, 4, 5),
                     I_1 = c(0, 0, 0, 0, 0, 0),
                     S_2 = c(0, 1, 2, 3, 4, 5),
                     I_2 = c(0, 0, 0, 0, 0, 0),
                     S_3 = c(0, 1, 2, 3, 4, 5),
                     I_3 = c(0, 0, 0, 0, 0, 0)),
                .Names = c("S_1", "I_1", "S_2", "I_2", "S_3", "I_3"),
                row.names = c(NA, -6L), class = "data.frame")

events <- structure(list(
    event      = c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3),
    time       = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    node       = c(2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6),
    dest       = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    n          = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
    proportion = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    select     = c(4, 5, 6, 4, 5, 6, 4, 5, 6, 4, 5, 6, 4, 5, 6),
    shift      = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)),
    .Names = c("event", "time", "node", "dest",
               "n", "proportion", "select", "shift"),
    row.names = c(NA, -15L), class = "data.frame")

model <- SISe3(u0        = u0,
               tspan     = 0:10,
               events    = events,
               phi       = rep(0, 6),
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

U <- structure(c(0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 2L,
                 0L, 2L, 0L, 2L, 0L, 3L, 0L, 3L, 0L, 3L, 0L, 4L, 0L, 4L, 0L, 4L,
                 0L, 5L, 0L, 5L, 0L, 5L, 0L, 15L, 0L, 15L, 0L, 15L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 15L, 0L, 15L,
                 0L, 15L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 15L, 0L, 15L, 0L, 15L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 15L, 0L, 15L, 0L, 15L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 15L, 0L, 15L,
                 0L, 15L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 15L, 0L, 15L, 0L, 15L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 15L, 0L, 15L, 0L, 15L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 15L, 0L, 15L,
                 0L, 15L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 15L, 0L, 15L, 0L, 15L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 15L, 0L, 15L, 0L, 15L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
               .Dim = c(36L, 11L))

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

## 6 Nodes
## 3 Age categories
## 2 Disease-states: Susceptible & Infected
##
## Zero probability of becoming infected.
##
## No individuals at t = 0
## At t = 1, all individuals enter in susceptible state
u0 <- structure(list(S_1 = c(0, 0, 0, 0, 0, 0),
                     I_1 = c(0, 0, 0, 0, 0, 0),
                     S_2 = c(0, 0, 0, 0, 0, 0),
                     I_2 = c(0, 0, 0, 0, 0, 0),
                     S_3 = c(0, 0, 0, 0, 0, 0),
                     I_3 = c(0, 0, 0, 0, 0, 0)),
                .Names = c("S_1", "I_1", "S_2", "I_2", "S_3", "I_3"),
                row.names = c(NA, -6L), class = "data.frame")

events <- structure(list(
    event      = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    time       = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    node       = c(2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6),
    dest       = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    n          = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
    proportion = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    select     = c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3),
    shift      = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)),
    .Names = c("event", "time", "node", "dest",
               "n", "proportion", "select", "shift"),
    row.names = c(NA, -15L), class = "data.frame")

model <- SISe3(u0        = u0,
               tspan     = 0:10,
               events    = events,
               phi       = rep(0, 6),
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

U <- structure(c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L,
                 0L, 1L, 0L, 2L, 0L, 2L, 0L, 2L, 0L, 3L, 0L, 3L, 0L, 3L, 0L, 4L,
                 0L, 4L, 0L, 4L, 0L, 5L, 0L, 5L, 0L, 5L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 1L, 0L, 1L, 0L, 1L, 0L, 2L, 0L, 2L, 0L, 2L, 0L, 3L, 0L, 3L,
                 0L, 3L, 0L, 4L, 0L, 4L, 0L, 4L, 0L, 5L, 0L, 5L, 0L, 5L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 2L, 0L, 2L, 0L, 2L,
                 0L, 3L, 0L, 3L, 0L, 3L, 0L, 4L, 0L, 4L, 0L, 4L, 0L, 5L, 0L, 5L,
                 0L, 5L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 2L,
                 0L, 2L, 0L, 2L, 0L, 3L, 0L, 3L, 0L, 3L, 0L, 4L, 0L, 4L, 0L, 4L,
                 0L, 5L, 0L, 5L, 0L, 5L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L,
                 0L, 1L, 0L, 2L, 0L, 2L, 0L, 2L, 0L, 3L, 0L, 3L, 0L, 3L, 0L, 4L,
                 0L, 4L, 0L, 4L, 0L, 5L, 0L, 5L, 0L, 5L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 1L, 0L, 1L, 0L, 1L, 0L, 2L, 0L, 2L, 0L, 2L, 0L, 3L, 0L, 3L,
                 0L, 3L, 0L, 4L, 0L, 4L, 0L, 4L, 0L, 5L, 0L, 5L, 0L, 5L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 2L, 0L, 2L, 0L, 2L,
                 0L, 3L, 0L, 3L, 0L, 3L, 0L, 4L, 0L, 4L, 0L, 4L, 0L, 5L, 0L, 5L,
                 0L, 5L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 2L,
                 0L, 2L, 0L, 2L, 0L, 3L, 0L, 3L, 0L, 3L, 0L, 4L, 0L, 4L, 0L, 4L,
                 0L, 5L, 0L, 5L, 0L, 5L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L,
                 0L, 1L, 0L, 2L, 0L, 2L, 0L, 2L, 0L, 3L, 0L, 3L, 0L, 3L, 0L, 4L,
                 0L, 4L, 0L, 4L, 0L, 5L, 0L, 5L, 0L, 5L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 1L, 0L, 1L, 0L, 1L, 0L, 2L, 0L, 2L, 0L, 2L, 0L, 3L, 0L, 3L,
                 0L, 3L, 0L, 4L, 0L, 4L, 0L, 4L, 0L, 5L, 0L, 5L, 0L, 5L, 0L),
               .Dim = c(36L, 11L))

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

## 6 Nodes
## 3 Age categories
## 2 Disease-states: Susceptible & Infected
##
## All individuals start in susceptible state, with a zero probability
## of becoming infected.
##
## At t = 3, all individuals exit.
u0 <- structure(list(S_1 = c(0, 1, 2, 3, 4, 5),
                     I_1 = c(0, 0, 0, 0, 0, 0),
                     S_2 = c(0, 1, 2, 3, 4, 5),
                     I_2 = c(0, 0, 0, 0, 0, 0),
                     S_3 = c(0, 1, 2, 3, 4, 5),
                     I_3 = c(0, 0, 0, 0, 0, 0)),
                .Names = c("S_1", "I_1", "S_2", "I_2", "S_3", "I_3"),
                row.names = c(NA, -6L), class = "data.frame")

events <- structure(list(
    event      = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    time       = c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3),
    node       = c(2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6),
    dest       = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    n          = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
    proportion = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    select     = c(4, 5, 6, 4, 5, 6, 4, 5, 6, 4, 5, 6, 4, 5, 6),
    shift      = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)),
    .Names = c("event", "time", "node", "dest",
               "n", "proportion", "select", "shift"),
    row.names = c(NA, -15L), class = "data.frame")

model <- SISe3(u0        = u0,
               tspan     = 0:10,
               events    = events,
               phi       = rep(0, 6),
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

U <- structure(c(0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 2L,
                 0L, 2L, 0L, 2L, 0L, 3L, 0L, 3L, 0L, 3L, 0L, 4L, 0L, 4L, 0L, 4L,
                 0L, 5L, 0L, 5L, 0L, 5L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L,
                 0L, 1L, 0L, 2L, 0L, 2L, 0L, 2L, 0L, 3L, 0L, 3L, 0L, 3L, 0L, 4L,
                 0L, 4L, 0L, 4L, 0L, 5L, 0L, 5L, 0L, 5L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 1L, 0L, 1L, 0L, 1L, 0L, 2L, 0L, 2L, 0L, 2L, 0L, 3L, 0L, 3L,
                 0L, 3L, 0L, 4L, 0L, 4L, 0L, 4L, 0L, 5L, 0L, 5L, 0L, 5L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
               .Dim = c(36L, 11L))

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

## 6 Nodes
## 3 Age categories
## 2 Disease-states: Susceptible & Infected
##
## All individuals start in susceptible state, with a zero probability of
## becoming infected.
##
## At t = 3, all individuals in age category 1 age.
## At t = 6, all individuals in age category 2 age.
u0 <- structure(list(S_1 = c(0, 1, 2, 3, 4, 5),
                     I_1 = c(0, 0, 0, 0, 0, 0),
                     S_2 = c(0, 1, 2, 3, 4, 5),
                     I_2 = c(0, 0, 0, 0, 0, 0),
                     S_3 = c(0, 1, 2, 3, 4, 5),
                     I_3 = c(0, 0, 0, 0, 0, 0)),
                .Names = c("S_1", "I_1", "S_2", "I_2", "S_3", "I_3"),
                row.names = c(NA, -6L), class = "data.frame")

events <- structure(list(event      = c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2),
                         time       = c(3, 3, 3, 3, 3, 6, 6, 6, 6, 6),
                         node       = c(2, 3, 4, 5, 6, 2, 3, 4, 5, 6),
                         dest       = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                         n          = c(1, 2, 3, 4, 5, 2, 4, 6, 8, 10),
                         proportion = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                         select     = c(4, 4, 4, 4, 4, 5, 5, 5, 5, 5),
                         shift      = c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2)),
                    .Names = c("event", "time", "node", "dest",
                        "n", "proportion", "select", "shift"),
                    row.names = c(NA, -10L), class = "data.frame")

model <- SISe3(u0        = u0,
               tspan     = 0:10,
               events    = events,
               phi       = rep(0, 6),
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

U <- structure(c(0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 2L,
                 0L, 2L, 0L, 2L, 0L, 3L, 0L, 3L, 0L, 3L, 0L, 4L, 0L, 4L, 0L, 4L,
                 0L, 5L, 0L, 5L, 0L, 5L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L,
                 0L, 1L, 0L, 2L, 0L, 2L, 0L, 2L, 0L, 3L, 0L, 3L, 0L, 3L, 0L, 4L,
                 0L, 4L, 0L, 4L, 0L, 5L, 0L, 5L, 0L, 5L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 1L, 0L, 1L, 0L, 1L, 0L, 2L, 0L, 2L, 0L, 2L, 0L, 3L, 0L, 3L,
                 0L, 3L, 0L, 4L, 0L, 4L, 0L, 4L, 0L, 5L, 0L, 5L, 0L, 5L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 2L, 0L, 1L, 0L, 0L, 0L, 4L, 0L, 2L,
                 0L, 0L, 0L, 6L, 0L, 3L, 0L, 0L, 0L, 8L, 0L, 4L, 0L, 0L, 0L, 10L,
                 0L, 5L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 2L, 0L, 1L, 0L, 0L,
                 0L, 4L, 0L, 2L, 0L, 0L, 0L, 6L, 0L, 3L, 0L, 0L, 0L, 8L, 0L, 4L,
                 0L, 0L, 0L, 10L, 0L, 5L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 2L, 0L, 1L, 0L, 0L, 0L, 4L, 0L, 2L, 0L, 0L, 0L, 6L, 0L, 3L, 0L,
                 0L, 0L, 8L, 0L, 4L, 0L, 0L, 0L, 10L, 0L, 5L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 3L, 0L, 0L, 0L, 0L, 0L, 6L, 0L, 0L,
                 0L, 0L, 0L, 9L, 0L, 0L, 0L, 0L, 0L, 12L, 0L, 0L, 0L, 0L, 0L,
                 15L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 3L, 0L, 0L,
                 0L, 0L, 0L, 6L, 0L, 0L, 0L, 0L, 0L, 9L, 0L, 0L, 0L, 0L, 0L, 12L,
                 0L, 0L, 0L, 0L, 0L, 15L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 3L, 0L, 0L, 0L, 0L, 0L, 6L, 0L, 0L, 0L, 0L, 0L, 9L, 0L,
                 0L, 0L, 0L, 0L, 12L, 0L, 0L, 0L, 0L, 0L, 15L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 3L, 0L, 0L, 0L, 0L, 0L, 6L, 0L, 0L,
                 0L, 0L, 0L, 9L, 0L, 0L, 0L, 0L, 0L, 12L, 0L, 0L, 0L, 0L, 0L,
                 15L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 3L, 0L, 0L,
                 0L, 0L, 0L, 6L, 0L, 0L, 0L, 0L, 0L, 9L, 0L, 0L, 0L, 0L, 0L, 12L,
                 0L, 0L, 0L, 0L, 0L, 15L, 0L),
               .Dim = c(36L, 11L))

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

## 6 Nodes
## 3 Age categories
## 2 Disease-states: Susceptible & Infected
##
## All individuals start in susceptible state, with a probability of
## becoming infected and then return to susceptible.
##
## No scheduled events
u0 <- structure(list(S_1 = c(0, 1, 2, 3, 4, 5),
                     I_1 = c(0, 0, 0, 0, 0, 0),
                     S_2 = c(0, 1, 2, 3, 4, 5),
                     I_2 = c(0, 0, 0, 0, 0, 0),
                     S_3 = c(0, 1, 2, 3, 4, 5),
                     I_3 = c(0, 0, 0, 0, 0, 0)),
                .Names = c("S_1", "I_1", "S_2", "I_2", "S_3", "I_3"),
                row.names = c(NA, -6L), class = "data.frame")

model <- SISe3(u0        = u0,
               tspan     = 0:10,
               events    = NULL,
               phi       = rep(1, 6),
               upsilon_1 = 1,
               upsilon_2 = 1,
               upsilon_3 = 1,
               gamma_1   = 1,
               gamma_2   = 1,
               gamma_3   = 1,
               alpha     = 1,
               beta_t1   = 1,
               beta_t2   = 1,
               beta_t3   = 1,
               beta_t4   = 1,
               end_t1    = 91,
               end_t2    = 182,
               end_t3    = 273,
               end_t4    = 365,
               epsilon   = 1)

result <- run(model, threads = 1, seed = 123L)
stopifnot(identical(model@G, result@G))
stopifnot(identical(model@S, result@S))
stopifnot(identical(sum(result@U[1:6,]), 0L))
stopifnot(all(apply(result@U[7:36,], 1, any)))
stopifnot(identical(model@ldata, result@ldata))
stopifnot(identical(model@sd, result@sd))
stopifnot(identical(model@tspan, result@tspan))
stopifnot(identical(model@u0, result@u0))
stopifnot(identical(model@events, result@events))

if (SimInf:::have_openmp()) {
    result_omp <- run(model, threads = 2, seed = 123L)
    stopifnot(identical(model@G, result_omp@G))
    stopifnot(identical(model@S, result_omp@S))
    stopifnot(identical(sum(result_omp@U[1:6,]), 0L))
    stopifnot(all(apply(result_omp@U[7:36,], 1, any)))
    stopifnot(identical(model@ldata, result_omp@ldata))
    stopifnot(identical(model@sd, result_omp@sd))
    stopifnot(identical(model@tspan, result_omp@tspan))
    stopifnot(identical(model@u0, result_omp@u0))
    stopifnot(identical(model@events, result_omp@events))
}

## Check extraction of seed parameter
.Call("SISe_run",
      demo_model(model = "SISe"),
      NULL,
      NULL,
      PACKAGE = "SimInf")

res <- .Call("SISe_run",
             demo_model(model = "SISe"),
             NULL,
             "1",
             PACKAGE = "SimInf")
stopifnot(identical(res$error, -6L))

.Call("SISe_run",
      demo_model(model = "SISe"),
      NULL,
      numeric(0),
      PACKAGE = "SimInf")

.Call("SISe_run",
      demo_model(model = "SISe"),
      NULL,
      integer(0),
      PACKAGE = "SimInf")

.Call("SISe_run",
      demo_model(model = "SISe"),
      NULL,
      1L,
      PACKAGE = "SimInf")

.Call("SISe_run",
      demo_model(model = "SISe"),
      NULL,
      1,
      PACKAGE = "SimInf")

res <- .Call("SISe_run",
             demo_model(model = "SISe"),
             NULL,
             NA_integer_,
             PACKAGE = "SimInf")
stopifnot(identical(res$error, -6L))

res <- .Call("SISe_run",
             demo_model(model = "SISe"),
             NULL,
             NA_real_,
             PACKAGE = "SimInf")
stopifnot(identical(res$error, -6L))

res <- .Call("SISe_run",
             demo_model(model = "SISe"),
             NULL,
             c(1L, 2L),
             PACKAGE = "SimInf")
stopifnot(identical(res$error, -6L))

res <- .Call("SISe_run",
             demo_model(model = "SISe"),
             NULL,
             c(1, 2),
             PACKAGE = "SimInf")
stopifnot(identical(res$error, -6L))

## Check extraction of number of threads
.Call("SISe_run",
      demo_model(model = "SISe"),
      NULL,
      NULL,
      PACKAGE = "SimInf")

.Call("SISe_run",
      demo_model(model = "SISe"),
      1L,
      NULL,
      PACKAGE = "SimInf")

.Call("SISe_run",
      demo_model(model = "SISe"),
      1,
      NULL,
      PACKAGE = "SimInf")

res <- .Call("SISe_run",
             demo_model(model = "SISe"),
             -1L,
             NULL,
             PACKAGE = "SimInf")
stopifnot(identical(res$error, -7L))

res <- .Call("SISe_run",
             demo_model(model = "SISe"),
             -1,
             NULL,
             PACKAGE = "SimInf")
stopifnot(identical(res$error, -7L))

res <- .Call("SISe_run",
             demo_model(model = "SISe"),
             "1",
             NULL,
             PACKAGE = "SimInf")
stopifnot(identical(res$error, -7L))

res <- .Call("SISe_run",
             demo_model(model = "SISe"),
             c(1L, 1L),
             NULL,
             PACKAGE = "SimInf")
stopifnot(identical(res$error, -7L))

res <- .Call("SISe_run",
             demo_model(model = "SISe"),
             c(1, 1),
             NULL,
             PACKAGE = "SimInf")
stopifnot(identical(res$error, -7L))

res <- .Call("SISe_run",
             demo_model(model = "SISe"),
             NA_integer_,
             NULL,
             PACKAGE = "SimInf")
stopifnot(identical(res$error, -7L))

res <- .Call("SISe_run",
             demo_model(model = "SISe"),
             NA_real_,
             NULL,
             PACKAGE = "SimInf")
stopifnot(identical(res$error, -7L))

## Check error codes in siminf_error
res <- tools::assertError(SimInf:::siminf_error("err"))
stopifnot(length(grep("'err' must be an integer vector of length 1",
                      res[[1]]$message)) > 0)

res <- tools::assertError(SimInf:::siminf_error(c(-1L, -1L)))
stopifnot(length(grep("'err' must be an integer vector of length 1",
                      res[[1]]$message)) > 0)

res <- tools::assertError(SimInf:::siminf_error(-1.1))
stopifnot(length(grep("'err' must be an integer vector of length 1",
                      res[[1]]$message)) > 0)

res <- tools::assertError(SimInf:::siminf_error(-10))
stopifnot(length(grep("Invalid model",
                      res[[1]]$message)) > 0)

res <- tools::assertError(SimInf:::siminf_error(-11))
stopifnot(length(grep("Unknown error code",
                      res[[1]]$message)) > 0)
