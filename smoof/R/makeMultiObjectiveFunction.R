#' Generator for multi-objective target functions.
#'
#' @template arg_name
#' @template arg_id
#' @template arg_description
#' @template arg_fn
#' @template arg_has_simple_signature
#' @template arg_par_set
#' @param n.objectives [\code{integer(1)}]\cr
#'   Number of objectives of the multi-objective function.
#' @template arg_noisy
#' @template arg_fn_mean
#' @param minimize [\code{logical}]\cr
#'   Logical vector of length \code{n.objectives} indicating if the corresponding
#'   objectives shall be minimized or maximized.
#'   Default is the vector with all components set to \code{TRUE}.
#' @template arg_vectorized
#' @template arg_constraint_fn
#' @param ref.point [\code{numeric}]\cr
#'   Optional reference point in the objective space, e.g., for hypervolume computation.
#' @return [\code{function}] Target function with additional stuff attached as attributes.
#' @examples
#' fn = makeMultiObjectiveFunction(
#'   name = "My test function",
#'   fn = function(x) c(sum(x^2), exp(x)),
#'   n.objectives = 2L,
#'   par.set = makeNumericParamSet("x", len = 1L, lower = -5L, upper = 5L)
#' )
#' print(fn)
#' @export
makeMultiObjectiveFunction = function(
  name = NULL,
  id = NULL,
  description = NULL,
  fn,
  has.simple.signature = TRUE,
  par.set,
  n.objectives,
  noisy = FALSE,
  fn.mean = NULL,
  minimize = rep(TRUE, n.objectives),
  vectorized = FALSE,
  constraint.fn = NULL,
  ref.point = NULL) {

  smoof.fn = makeObjectiveFunction(
    name, id, description, fn,
    has.simple.signature, par.set, n.objectives,
    noisy, fn.mean, minimize, vectorized, constraint.fn
  )

  if (!is.null(ref.point)) {
    assertNumeric(ref.point, len = n.objectives, any.missing = FALSE, all.missing = FALSE)
  }

  smoof.fn = setAttribute(smoof.fn, "ref.point", ref.point)
  class(smoof.fn) = c("smoof_multi_objective_function", class(smoof.fn))

  return(smoof.fn)
}
