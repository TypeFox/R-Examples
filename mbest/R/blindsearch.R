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


# Perform a "blind" linesearch: don't check any curvature or
# sufficient decrease conditions; always declare the search to have
# converged after one iteration.
#


# default search parameters
blindsearch.control <- function()
{
    list()
}


# Start a line search, given initial value, derivative and step size.
# Return value includes fields `step` giving next step to true, and
# `converged` indicating whether or not search has converged.
blindsearch <- function(value, deriv, step, control = list())
{
    control <- do.call("blindsearch.control", control)

    if (!is.numeric(value) || is.na(value))
        stop("missing initial value")
    if (!is.numeric(deriv) || is.na(deriv))
        stop("missing initial derivative")
    if (!(deriv < 0 && is.finite(deriv)))
        stop("initial derivative must be negative and finite")

    best <- NULL

    z <- list(step = step, value0 = value, deriv0 = deriv, 
              best = best, converged = FALSE, control=control)
    class(z) <- "blindsearch"
    z
}


print.blindsearch <- function(x, digits = max(3L, getOption("digits") - 3L),
                              ...)
{
    if (x$converged) {
        cat("Converged blind linesearch (step = ",
            format(x$step, digits), ")\n", sep="")
    } else {
        cat("Blind search in progress; next step = ",
            format(x$step, digits), "\n", sep="")
    }
}


# update search with new value and derivative
update.blindsearch <- function(object, value, deriv, ...)
{
    step <- object$step

    converged <- FALSE

    if (is.na(value))
        stop("function value is NA")
    if (is.na(deriv))
        stop("function derivative is NA")

    test <- list(position = object$step, value = value, deriv = deriv)

    object$best <- test
    object$converged <- TRUE
    object$value <- value
    object$deriv <- deriv

    object
}
