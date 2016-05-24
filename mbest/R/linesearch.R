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


# Perform a one-dimensional line search to find a step length satisfying
# the strong Wolfe conditions (sufficient decrease and curvature condition).
# Useful as part of a function minimization algorithm.  
#
# This module implements the algorithm described by More and Thuente (1994),
# which is guaranteed to converge after a finite number of iterations.
#
# References:
#
# * More, J. J. and Thuente, D. J. (1994). Line search algorithms
#   with guaranteed sufficient decrease. /ACM Transactions on Mathematical
#   Software/ 20(3):286-307.
#   <http://doi.acm.org/10.1145/192115.192132>
#
# * Nocedal, J. and Wright, S. J. (2006). /Numerical Optimization/, 2nd ed.
#   Springer.
#   <http://www.springer.com/mathematics/book/978-0-387-30303-1>
#


# default search parameters
linesearch.control <- function(value.tol = 1e-4, deriv.tol = 0.9,
                               step.tol = 1e-7, step.min = 1e-10,
                               step.max = 1e10, extrap.lower = 1.1,
                               extrap.upper = 4.0, extrap.max = 0.66,
                               bisection.width = 0.66)
{
    if (!is.numeric(value.tol) || value.tol <= 0)
        stop("value of 'value.tol' must be > 0")
    if (!is.numeric(deriv.tol) || deriv.tol <= 0)
        stop("value of 'deriv.tol' must be > 0")
    if (!is.numeric(step.tol) || step.tol <= 0)
        stop("value of 'step.tol' must be > 0")
    if (!is.numeric(step.min) || step.min <= 0)
        stop("value of 'step.min' must be > 0")
    if (!is.numeric(step.max) || step.max <= step.min)
        stop("value of 'step.max' must be > step.min")
    if (!is.numeric(extrap.lower) || !(extrap.lower > 1))
        stop("value of 'extrap.lower' must be > 1")
    if (!is.numeric(extrap.upper) || !(extrap.upper > extrap.lower))
        stop("value of 'extrap.upper' must be > extrap.lower")
    if (!is.numeric(extrap.max) || !(0 < extrap.max && extrap.max < 1))
        stop("value of 'extrap.max' must be in range (0,1)")
    if (!is.numeric(bisection.width)
            || !(0 < bisection.width && bisection.width < 1))
        stop("value of 'bisection.width' must be in range (0,1)")

    list(value.tol = value.tol, deriv.tol = deriv.tol,
         step.tol = step.tol, step.min = step.min, step.max = step.max,
         extrap.lower = extrap.lower, extrap.upper = extrap.upper,
         extrap.max = extrap.max, bisection.width = bisection.width)
}


# Start a line search, given initial value, derivative and step size.
# Return value includes fields `step` giving next step to true, and
# `converged` indicating whether or not search has converged.
linesearch <- function(value, deriv, step, control = list())
{
    control <- do.call("linesearch.control", control)

    if (!is.numeric(value) || is.na(value))
        stop("missing initial value")
    if (!is.numeric(deriv) || is.na(deriv))
        stop("missing initial derivative")
    if (!(deriv < 0 && is.finite(deriv)))
        stop("initial derivative must be negative and finite")
    if (!(control$step.min < step && step < control$step.max))
        stop("initial step must be in range (step.min, step.max)")

    gtest <- deriv * control$value.tol
    ftest <- value + step * gtest
    lower <- list(position = 0, value = value, deriv = deriv)
    upper <- lower
    interval <- c(0, step + control$extrap.upper * step)

    z <- list(step = step, value0 = value, deriv0 = deriv, ftest = ftest,
              gtest = gtest, bracket = FALSE, auxfun = TRUE, lower = lower,
              upper = upper, interval = interval, old.width1 = Inf,
              old.width2 = Inf, converged = FALSE, control=control)
    class(z) <- "linesearch"
    z
}


print.linesearch <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
    if (x$converged) {
        cat("Converged linesearch (step = ",
            format(x$step, digits), ")\n", sep="")
    } else {
        cat("Linesearch in progress; next step = ",
            format(x$step, digits), "\n", sep="")
    }
}


# update search with new value and derivative
update.linesearch <- function(object, value, deriv, ...)
{
    control <- object$control
    interval <- object$interval
    ftest <- object$ftest
    gtest <- object$gtest
    step <- object$step

    converged <- FALSE
    stuck <- FALSE

    if (is.na(value))
        stop("function value is NA")
    if (is.na(deriv))
        stop("function derivative is NA")

    if (value <= ftest && abs(deriv) <= -(control$deriv.tol * object$deriv0)) {
        converged <- TRUE
    } else if (object$bracket && step <= interval[1] || step >= interval[2]) {
        warning("lack of numerical precision prevents further progress")
        stuck <- TRUE
    } else if (object$bracket
               && diff(interval) <= control$step.tol * interval[1]) {
        warning("step size within tolerance")
        stuck <- TRUE
    } else if (step == control$step.max && value <= ftest && deriv <= gtest) {
        warning("at step.max")
        stuck <- TRUE
    } else if (step == control$step.min && (value > ftest || deriv >= gtest)) {
        warning("at step.min")
        stuck <- TRUE
    }

    if (converged || stuck) {
        object$converged <- converged
        object$value <- value
        object$deriv <- deriv
        return(object)
    }

    test <- list(position = object$step, value = value, deriv = deriv)
    update.linesearch.step(object, test)
}


update.linesearch.step <- function(object, test)
{
    ftest <- object$ftest
    gtest <- object$gtest

    # Use auxiliary function to start with (psi, p.290)
    # If the auxiliary function is nonpositive and the derivative of the
    # original function is positive, switch to the original function (p.298)
    aux <- object$auxfun && !(test$value <= ftest && test$deriv > 0)

    test0 <- test
    lower0 <- object$lower
    upper0 <- object$upper

    if (aux) {
        modify <- function(eval) {
            list(position = eval$position,
                 value = (eval$value) - (eval$position) * gtest,
                 deriv = (eval$deriv) - gtest)
        }

        test <- modify(test0)
        object$lower <- modify(lower0)
        object$upper <- modify(upper0)
    } else {
        object$auxfun <- FALSE
    }

    object <- unsafe.step(object, test)

    if (aux) {
        restore <- function(eval) {
            if (eval$position == test0$position) {
                test0
            } else if (eval$position == lower0$position) {
                lower0
            } else {
                upper0
            }
        }

        object$lower <- restore(object$lower)
        object$upper <- restore(object$upper)
    }

    maybe.bisect(object)
}



# If the length of the interval doesn't decrease by delta after two
# iterations, then perform bisection instead (pp.292-293).
maybe.bisect <- function(object)
{

    if (!object$bracket)
        return(object)

    l <- object$lower$position
    u <- object$upper$position
    w <- diff(object$interval)
    control <- object$control

    if (w >= (control$bisection.width) * (object$old.width2)) {
        object$step <- l + 0.5 * (u - l)
        object$ftest <- object$value0 + object$step * object$gtest
    }

    object$old.width2 <- object$old.width1
    object$old.width1 <- w

    object
}


unsafe.step <- function(object, test)
{
    lower <- object$lower
    upper <- object$upper
    gtest <- object$gtest
    control <- object$control
    bracket <- object$bracket

    # compute new trial step
    step <- trial.value(object$interval, lower, upper, test, bracket, control)

    # update the search interval
    upd <- update.interval(lower, upper, test, step, bracket, control)

    # safeguard the step
    step <- max(control$step.min, min(control$step.max, step))
    ftest <- object$value0 + step * gtest

    object$bracket <- upd$bracket
    object$lower <- upd$lower
    object$upper <- upd$upper
    object$interval <- upd$interval
    object$step <- step
    object$ftest <- ftest

    object
}


# (Modified) Updating Algorithm (pp. 291, 297-298)
update.interval <- function(lower, upper, test, step, bracket, control)
{
    fl <- lower$value
    gl <- lower$deriv
    t <- test$position
    ft <- test$value
    gt <- test$deriv

    t1 <- step

    bracket.endpoints <- function(lower, upper) {
        list(bracket = TRUE,
             interval = sort(c(lower$position, upper$position)),
             lower = lower,
             upper = upper)
    }

    # Case U1: higher function value
    if (ft > fl) {
        bracket.endpoints(lower, test)

    # Case U2: lower function value, derivatives same sign
    } else if (sign(gt) == sign(gl)) {
        if (bracket) {
            bracket.endpoints(test, upper)

        # If we haven't found a suitable interval yet, then we
        # extrapolate (p.291)
        } else {
            list(bracket = FALSE,
                 interval = c(t1 + (control$extrap.lower) * (t1 - t),
                              t1 + (control$extrap.upper) * (t1 - t)),
                 lower = test,
                 upper = upper)
        }

    # Case U3: lower function value, derivatives same sign
    } else {
        bracket.endpoints(test, lower)
    }
}


# Trial value selection (Sec. 4, pp. 298-300)
trial.value <- function(interval, lower, upper, test, bracket, control)
{
    tmin <- interval[1]
    tmax <- interval[2]

    l <- lower$position
    fl <- lower$value
    gl <- lower$deriv

    u <- upper$position
    fu <- upper$value
    gu <- upper$deriv

    t <- test$position
    ft <- test$value
    gt <- test$deriv

    c <- suppressWarnings(cubic.min(l, fl, gl, t, ft, gt)) # may be NaN
    q <- quadr.min(l, fl, gl, t, ft)
    s <- secant.min(l, gl, t, gt)


    # Case 1: higher function value
    if (ft > fl) {
        if (abs(c - l) < abs (q - l)) {
            c
        } else {
            c + 0.5 * (q - c)
        }

    # Case 2: lower function value, derivative opposite sign
    } else if (sign(gt) != sign(gl)) {
        if (abs(c - t) >= abs(s - t)) {
            c
        } else {
            s
        }

    # Case 3: lower function value, derivatives same sign, lower derivative
    } else if (abs(gt) <= abs(gl)) {
        # if the cubic tends to infinity in the direction of the
        # step and the minimum of the cubic is beyond t...
        c1 <- if (!is.nan(c) && sign(c - t) == sign(c - l)) {
                 # ...then the cubic step is ok
                 c
              # ...otherwise replace it with the endpoint in the
              # direction of t (secant step in paper)
              } else if (t > l) tmax else tmin

        if (bracket) {
            guard.extrap <- if (t > l)
                 function(x) min(x, t + control$extrap.max * (u - t))
            else function(x) max(x, t + control$extrap.max * (u - t))

           if (abs(t - c1) < abs(t - s))
               guard.extrap(c1)
           else
               guard.extrap(s)
        } else {
            clip <- function(x) max(tmin, min(tmax, x))

            #extrapolate to farthest of cubic and secant steps
            if (abs(t - c1) > abs(t - s)) clip(c1)
            else clip(s)
        }

    # Case 4: lower function value, derivatives same sign, higher derivative
    } else {
        if (bracket) {
            cubic.min(t, ft, gt, u, fu, gu)
        } else if (t > l) {
            tmax
        } else {
            tmin
        }
    }
}


quadr.min <- function(u, fu, du, v, fv)
{
    a <- v - u
    u + (du / ((fu - fv) / a + du) / 2) * a
}


secant.min <- function(u, du, v, dv)
{
    a <- u - v
    v + dv / (dv - du) * a
}


cubic.min <- function(u, fu, du, v, fv, dv)
{
    d <- v - u
    theta <- (fu - fv) * 3 / d + du + dv
    s <- max(abs(c(theta, du, dv)))
    a <- theta / s
    gamma0 <- s * sqrt(a * a - (du / s) * (dv / s))
    gamma <- if (v < u) -gamma0 else gamma0
    p <- gamma - du + theta
    q <- gamma - du + gamma + dv
    r <- p / q
    u + r * d
}
