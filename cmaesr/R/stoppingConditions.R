#' @title Stopping condition: maximal iterations.
#'
#' @description Stop on maximal number of iterations.
#'
#' @param max.iter [integer(1)]\cr
#'   Maximal number of iterations.
#'   Default is \code{100}.
#' @return [\code{cma_stopping_condition}]
#' @family stopping conditions
#' @export
stopOnMaxIters = function(max.iter = 100L) {
  assertInt(max.iter, na.ok = FALSE)
  force(max.iter)
  return(makeStoppingCondition(
    name = "maxIter",
    message = sprintf("MaxIter: reached maximal number of iterations/generations %i.", max.iter),
    stop.fun = function(envir = parent.frame()) {
      return(envir$iter > max.iter)
    }
  ))
}

# @title Stopping condition: indefinite covariance matrix.
#
# @description Stop if covariance matrix is not positive definite anymore.
#
# @return [\code{cma_stopping_condition}]
# @family stopping conditions
#NOTE: this one is not exported. However, it is always prepended to the list of
# stopping conditions, since
stopOnIndefCovMat = function() {
  return(makeStoppingCondition(
    name = "indefCovMat",
    message = "Covariance matrix is not numerically positive definite.",
    stop.fun = function(envir = parent.frame()) {
      e.values = envir$e$values
      return(any(is.na(e.values)) || any(e.values <= sqrt(.Machine$double.eps) * abs(e.values[1L])))
    }
  ))
}

#' @title Stopping condition: optimal params.
#'
#' @description Stop if euclidean distance of parameter is below
#' some tolerance value.
#'
#' @param opt.param [\code{numeric}]\cr
#'   Known optimal parameter settings.
#' @param tol [\code{numeric(1)}]\cr
#'   Tolerance value.
#'   Default is \eqn{1e^{-8}}.
#' @return [\code{cma_stopping_condition}]
#' @family stopping conditions
#' @export
stopOnOptParam = function(opt.param, tol = 1e-8) {
  assertNumeric(opt.param, any.missing = FALSE, all.missing = FALSE)
  assertNumber(tol, lower = 0, na.ok = FALSE, finite = TRUE)
  force(opt.param)
  force(tol)
  return(makeStoppingCondition(
    name = "optParamTol",
    message = sprintf("Optimal parameters approximated nicely (gap < %.2f).", tol),
    stop.fun = function(envir = parent.frame()) {
      return(sqrt(sum(envir$best.param - opt.param)^2) < tol)
    }
  ))
}

#' @title Stopping condition: optimal objective value.
#'
#' @description Stop if best solution is close to optimal objective value.
#'
#' @param opt.value [\code{numeric(1)}]\cr
#'   Known optimal objective function value.
#' @param tol [\code{numeric(1)}]\cr
#'   Tolerance value.
#'   Default is \eqn{1e^{-8}}.
#' @return [\code{cma_stopping_condition}]
#' @family stopping conditions
#' @export
stopOnOptValue = function(opt.value, tol = 1e-8) {
  assertNumber(opt.value, na.ok = FALSE)
  assertNumber(tol, lower = 0, na.ok = FALSE, finite = TRUE)
  force(opt.value)
  force(tol)
  return(makeStoppingCondition(
    name = "optValTol",
    message = sprintf("Optimal function value approximated nicely (gap < %.10f).", tol),
    stop.fun = function(envir = parent.frame()) {
      return(abs(envir$best.fitness - opt.value) < tol)
    }
  ))
}

#' @title Stopping condition: maximal time.
#'
#' @description Stop if maximal running time budget is reached.
#'
#' @param budget [\code{integer(1)}]\cr
#'   Time budget in seconds.
#' @return [\code{cma_stopping_condition}]
#' @family stopping conditions
#' @export
stopOnTimeBudget = function(budget) {
  assertInt(budget, na.ok = FALSE, lower = 1L)
  force(budget)
  return(makeStoppingCondition(
    name = "timeBudget",
    message = sprintf("Time budget of %i [secs] reached.", budget),
    stop.fun = function(envir = parent.frame()) {
      return(difftime(Sys.time(), envir$start.time, units = "secs") > budget)
    }
  ))
}

#' @title Stopping condition: maximal funtion evaluations.
#'
#' @description Stop if maximal number of function evaluations is reached.
#'
#' @param max.evals [\code{integer(1)}]\cr
#'   Maximal number of allowed function evaluations.
#' @return [\code{cma_stopping_condition}]
#' @export
stopOnMaxEvals = function(max.evals) {
  assertInt(max.evals, na.ok = FALSE, lower = 1L)
  force(max.evals)
  return(makeStoppingCondition(
    name = "maxEvals",
    message = sprintf("Maximal number of %i function evaluations reached.", max.evals),
    stop.fun = function(envir = parent.frame()) {
      return(envir$n.evals >= max.evals)
    }
  ))
}

#' @title Stopping condition: low standard deviation.
#'
#' @description Stop if the standard deviation falls below a tolerance value
#' in all coordinates?
#'
#' @param tol [\code{integer(1)}]\cr
#'   Tolerance value.
#' @return [\code{cma_stopping_condition}]
#' @export
#FIXME: default value is 1e-12 * sigma. Here we have no access to the sigma value.
stopOnTolX = function(tol = 1e-12) {
  assertInt(tol, na.ok = FALSE)
  force(tol)
  return(makeStoppingCondition(
    name = "tolX",
    message = sprintf("Standard deviation below tolerance in all coordinates."),
    stop.fun = function(envir = parent.frame()) {
      return(all(envir$D < tol) && all((envir$sigma * envir$p.c) < tol))
    }
  ))
}

#' @title Stopping condition: principal axis.
#'
#' @description Stop if addition of 0.1 * sigma in a principal axis
#' direction does not change mean value.
#'
#' @return [\code{cma_stopping_condition}]
#' @family stopping conditions
#' @export
stopOnNoEffectAxis = function() {
  return(makeStoppingCondition(
    name = "noEffectAxis",
    message = "Addition of 0.1 times sigma does not change mean value.",
    stop.fun = function(envir = parent.frame()) {
      ii = (envir$iter %% envir$n) + 1L
      ui = envir$e$vectors[, ii]
      lambdai = sqrt(envir$e$values[ii])
      m = envir$m
      return(sum((m - (m + 0.1 * envir$sigma * lambdai * ui))^2) < .Machine$double.eps)
    }
  ))
}

#' @title Stopping condition: standard deviation in coordinates.
#'
#' @description Stop if addition of 0.2 * standard deviations in any
#' coordinate does not change mean value.
#'
#' @return [\code{cma_stopping_condition}]
#' @family stopping conditions
#' @export
stopOnNoEffectCoord = function() {
  return(makeStoppingCondition(
    name = "noEffectCoord",
    message = "Addition of 0.2 times sigma in any coordinate does not change mean value.",
    stop.fun = function(envir = parent.frame()) {
      m = envir$m
      return(sum((m - (m + 0.2 * envir$sigma))^2) < .Machine$double.eps)
    }
  ))
}

#' @title Stopping condition: high condition number.
#'
#' @description Stop if condition number of covariance matrix exceeds
#' tolerance value.
#'
#' @param tol [\code{numeric(1)}]\cr
#'   Tolerance value.
#'   Default is \code{1e14}.
#' @return [\code{cma_stopping_condition}]
#' @family stopping conditions
#' @export
stopOnCondCov = function(tol = 1e14) {
  assertNumber(tol, na.ok = FALSE, lower = 0, finite = TRUE)
  force(tol)
  return(makeStoppingCondition(
    name = "conditionCov",
    message = sprintf("Condition number of covariance matrix exceeds %f", tol),
    stop.fun = function(envir = parent.frame()) {
      return(kappa(envir$C) > tol)
    }
  ))
}
