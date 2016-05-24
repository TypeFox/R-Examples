#' Compute the Expected Running Time (ERT) performance measure.
#'
#' The functions can be called in two different ways
#' \itemize{
#'   \item{1. Pass a vector of function evaluations and a logical vector which
#'   indicates which runs were successful (see details).}
#'   \item{2. Pass a vector of function evaluation, a vector of reached target
#'   values and a single target value. In this case the logical vector of
#'   option 1. is computed internally.}
#' }
#'
#' @details The Expected Running Time (ERT) is one of the most popular performance
#' measures in optimization. It is defined as the expected number of function
#' evaluations needed to reach a given precision level, i. e., to reach a certain
#' objective value for the first time.
#'
#' @references A. Auger and N. Hansen. Performance evaluation of an advanced local
#' search evolutionary algorithm. In Proceedings of the IEEE Congress on Evolutionary
#' Computation (CEC 2005), pages 1777-1784, 2005.
#'
#' @param fun.evals [\code{numeric}]\cr
#'   Vector containing the number of function evaluations.
#' @param fun.success.runs [\code{logical}]\cr
#'   Boolean vector indicating which algorithm runs were successful,
#'   i. e., which runs reached the desired target value. Default is \code{NULL}.
#' @param fun.reached.target.values [\code{numeric} | \code{NULL}]\cr
#'   Numeric vector with the objective values reached in the runs. Default is
#'   \code{NULL}.
#' @param fun.target.value [\code{numeric(1)} | \code{NULL}]\cr
#'   Target value which shall be reached. Default is \code{NULL}.
#' @param penalty.value [\code{numeric(1)}]\cr
#'   Penalty value which should be returned if none of the algorithm runs
#'   was successful. Default is \code{Inf}.
#' @return [\code{numeric(1)}]
#'   Estimated Expected Running Time.
#' @export
computeExpectedRunningTime = function(fun.evals,
  fun.success.runs = NULL,
  fun.reached.target.values = NULL,
  fun.target.value = NULL,
  penalty.value = Inf) {
  #FIXME: maybe enable missing values and offer inpute mechanism?
  assertInteger(fun.evals, lower = 1L, any.missing = FALSE)
  assertNumber(penalty.value, na.ok = FALSE)
  n = length(fun.evals)

  # sanity check that one of the options is used (see docs).
  if (!xor(!is.null(fun.success.runs), (!is.null(fun.reached.target.values) || !is.null(fun.target.value)))) {
    stopf("Either 'fun.success.runs' or 'fun.reached.target.values' and 'fun.target.value' need to be specified,
      but not both or none.")
  }

  # compute successful runs
  if (!is.null(fun.reached.target.values)) {
    if (is.null(fun.target.value)) {
      stopf("You need to pass a 'fun.target.value' in case you passed 'fun.reached.target.values'.")
    }
    assertNumeric(fun.reached.target.values, len = n, any.missing = FALSE)
    assertNumber(fun.target.value, na.ok = FALSE)
    fun.success.runs = (fun.reached.target.values <= fun.target.value)
  }

  assertLogical(fun.success.runs, len = n, any.missing = FALSE)

  # Finally compute the ERT
  # compute the success rate (since fun.success.runs is logical, this is correct)
  ps = mean(fun.success.runs)

  # average number of function evaluations for successful runs
  RT.succ = mean(fun.evals[fun.success.runs])

  # average number of function evaluations for unsuccessful runs
  RT.fail = mean(fun.evals[!fun.success.runs])

  if (ps == 1.0) {
    return(RT.succ)
  } else if (ps == 0.0) {
    return(penalty.value)
  } else {
    return(RT.succ + (1 - ps) / ps * RT.fail)
  }
}
