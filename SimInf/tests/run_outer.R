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

## Define some variables
tol = 1e-8
x = seq(from = 0.95, to = 1.05, by = 0.01)
y = seq(from = 0.95, to = 1.05, by = 0.01)

## Check gdata names
model <- demo_model()
names(model@gdata) <- NULL
res <- tools::assertError(
    run_outer(x, y, model, alpha ~ upsilon, function(model) 1))
stopifnot(length(grep("'names[(]model@gdata[)]' is NULL",
                      res[[1]]$message)) > 0)

## Check formula argument
res <- tools::assertError(run_outer(x, y, demo_model(), NULL, function(model) 1))
stopifnot(length(grep("'formula' argument is NULL",
                      res[[1]]$message)) > 0)

## Check FUN argument
res <- tools::assertError(run_outer(x, y, demo_model(), a ~ b, NULL))
stopifnot(length(grep("'FUN' argument is NULL",
                      res[[1]]$message)) > 0)

## Check lhs
res <- tools::assertError(
    run_outer(x, y, demo_model(),  ~ upsilon, function(model) 1))
stopifnot(length(grep("Invalid parameters on the left side of the formula",
                      res[[1]]$message)) > 0)

res <- tools::assertError(
    run_outer(x, y, demo_model(), alph ~ upsilon, function(model) 1))
stopifnot(length(grep("Unmatched parameters on the left hand side of the formula",
                      res[[1]]$message)) > 0)

## Check rhs
res <- tools::assertError(
    run_outer(x, y, demo_model(), alpha ~ upsilon:alpha, function(model) 1))
stopifnot(length(grep("Invalid parameters on the right side of the formula",
                      res[[1]]$message)) > 0)

res <- tools::assertError(
    run_outer(x, y, demo_model(), alph ~ upsilon, function(model) 1))
stopifnot(length(grep("Unmatched parameters on the left hand side of the formula",
                      res[[1]]$message)) > 0)

## Check run_outer
z_exp <- structure(
    c(0.008745225, 0.00883728, 0.008929335, 0.00902139,
      0.009113445, 0.0092055, 0.009297555, 0.00938961, 0.009481665,
      0.00957372, 0.009665775, 0.00883728, 0.008930304, 0.009023328,
      0.009116352, 0.009209376, 0.0093024, 0.009395424, 0.009488448,
      0.009581472, 0.009674496, 0.00976752, 0.008929335, 0.009023328,
      0.009117321, 0.009211314, 0.009305307, 0.0093993, 0.009493293,
      0.009587286, 0.009681279, 0.009775272, 0.009869265, 0.00902139,
      0.009116352, 0.009211314, 0.009306276, 0.009401238, 0.0094962,
      0.009591162, 0.009686124, 0.009781086, 0.009876048, 0.00997101,
      0.009113445, 0.009209376, 0.009305307, 0.009401238, 0.009497169,
      0.0095931, 0.009689031, 0.009784962, 0.009880893, 0.009976824,
      0.010072755, 0.0092055, 0.0093024, 0.0093993, 0.0094962, 0.0095931,
      0.00969, 0.0097869, 0.0098838, 0.0099807, 0.0100776, 0.0101745,
      0.009297555, 0.009395424, 0.009493293, 0.009591162, 0.009689031,
      0.0097869, 0.009884769, 0.009982638, 0.010080507, 0.010178376,
      0.010276245, 0.00938961, 0.009488448, 0.009587286, 0.009686124,
      0.009784962, 0.0098838, 0.009982638, 0.010081476, 0.010180314,
      0.010279152, 0.01037799, 0.009481665, 0.009581472, 0.009681279,
      0.009781086, 0.009880893, 0.0099807, 0.010080507, 0.010180314,
      0.010280121, 0.010379928, 0.010479735, 0.00957372, 0.009674496,
      0.009775272, 0.009876048, 0.009976824, 0.0100776, 0.010178376,
      0.010279152, 0.010379928, 0.010480704, 0.01058148, 0.009665775,
      0.00976752, 0.009869265, 0.00997101, 0.010072755, 0.0101745,
      0.010276245, 0.01037799, 0.010479735, 0.01058148, 0.010683225
      ), .Dim = c(11L, 11L))

x = seq(from = 0.95, to = 1.05, by = 0.01)
y = seq(from = 0.95, to = 1.05, by = 0.01)
run_f <- function(model, N) {
    model@gdata["upsilon"] * model@gdata["beta_t1"] * N
}
z_obs <- run_outer(x, y, demo_model(), upsilon ~ beta_t1, run_f, N = 3)

stopifnot(all(abs(z_obs - z_exp) < tol))
