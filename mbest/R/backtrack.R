# Copyright 2014 Patrick O. Perry
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


# Perform a one-dimensional backtracking line search to find a step length
# satisfying the Armijo condition (sufficient decrease).
#
# References:
#
# * Nocedal, J. and Wright, S. J. (2006). /Numerical Optimization/, 2nd ed.
#   Springer.
#   <http://www.springer.com/mathematics/book/978-0-387-30303-1>
#


# default search parameters
backtrack.control <- function(value.tol = 1e-4, contract = 0.5,
                              step.min = 1e-10)
{
    if (!is.numeric(value.tol) || value.tol <= 0)
        stop("value of 'value.tol' must be > 0")
    if (!is.numeric(contract) || !(0 < contract && contract < 1))
        stop("value  of 'contract' must be in (0,1)")
    if (!is.numeric(step.min) || step.min <= 0)
        stop("value of 'step.min' must be > 0")

    list(value.tol = value.tol, contract = contract, step.min = step.min)
}


# Start a line search, given initial value, derivative and step size.
# Return value includes fields `step` giving next step to true, and
# `converged` indicating whether or not search has converged.
backtrack <- function(value, deriv, step, control = list())
{
    control <- do.call("backtrack.control", control)

    if (!is.numeric(value) || is.na(value))
        stop("missing initial value")
    if (!is.numeric(deriv) || is.na(deriv))
        stop("missing initial derivative")
    if (!(deriv < 0 && is.finite(deriv)))
        stop("initial derivative must be negative and finite")
    if (!(control$step.min < step))
        stop("initial step must be > step.min")

    gtest <- deriv * control$value.tol
    ftest <- value + step * gtest
    best <- list(position = 0, value = value, deriv = deriv)

    z <- list(step = step, value0 = value, deriv0 = deriv, ftest = ftest,
              gtest = gtest, best = best, converged = FALSE, control=control)
    class(z) <- "backtrack"
    z
}


print.backtrack <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
    if (x$converged) {
        cat("Converged backtracking linesearch (step = ",
            format(x$step, digits), ")\n", sep="")
    } else {
        cat("Backtracking linesearch in progress; next step = ",
            format(x$step, digits), "\n", sep="")
    }
}


# update search with new value and derivative
update.backtrack <- function(object, value, deriv, ...)
{
    control <- object$control
    ftest <- object$ftest
    gtest <- object$gtest
    step <- object$step

    converged <- FALSE
    stuck <- FALSE

    if (is.na(value))
        stop("function value is NA")
    if (is.na(deriv))
        stop("function derivative is NA")

    test <- list(position = object$step, value = value, deriv = deriv)
    if (object$best$value > value)
        object$best <- test

    if (value <= ftest) {
        converged <- TRUE
    } else if (step == control$step.min) {
        warning("at step.min")
        stuck <- TRUE
    }

    if (converged || stuck) {
        object$converged <- converged
        object$value <- value
        object$deriv <- deriv
        return(object)
    }

    update.backtrack.step(object, test)
}


update.backtrack.step <- function(object, test)
{
    control <- object$control

    # compute new trial step
    step <- (test$position) * (control$contract)

    # safeguard the step
    step <- max(control$step.min, step)

    object$ftest <- object$value0 + step * (object$gtest)
    object$step <- step

    object
}
