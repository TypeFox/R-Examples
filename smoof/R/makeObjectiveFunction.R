# Internal generator for function of smoof type.
#
# @param name [\code{character(1)}]\cr
#   Optional function name used e.g. in plots.
# @param description [\code{character(1)} | \code{NULL}]\cr
#   Optional function description.
# @param fn [\code{function}]\cr
#   Target function.
# @param has.simple.signature [\code{logical(1)}]\cr
#   Set this to \code{TRUE} if the target function expects a vector as input and \code{FALSE}
#   if it expects a named list of values. The latter is needed if the function depends on mixed
#   parameters. Default is \code{TRUE}.
# @param par.set [\code{\link[ParamHelpers]{ParamSet}}]\cr
#   Parameter set describing different ascpects of the target function parameters, i. e.,
#   names, lower and/or upper bounds, types and so on. See \code{\link[ParamHelpers]{makeParamSet}}
#   for further information.
# @param n.objectives [\code{integer(1)}]\cr
#   Number of objectives of the multi-objective function.
# @param noisy [\code{logical(1)}]\cr
#   Is the function noisy? Defaults to \code{FALSE}.
# @param fn.mean [\code{function}]\cr
#   Optional true mean function in case of a noisy objective function. This functions should
#   have the same mean as \code{fn}.
# @param minimize [\code{logical}]\cr
#   Logical vector of length \code{n.objectives} indicating which objectives shall
#   be minimzed/maximized.
#   The default is \code{TRUE} \code{n.objectives} times.
# @param vectorized [\code{logical(1)}]\cr
#   Can the handle \dQuote{vector} input, i. e., does it accept matrix of
#   parameters? Default is \code{FALSE}.
# @param constraint.fn [\code{function | NULL}]\cr
#   Function which returns a logical vector indicating which indicates whether certain conditions
#   are met or not. Default is \code{NULL}, which means, that there are no constraints (beside possible)
#   box constraints defined via the \code{par.set} argument.
# @return [\code{function}] Target function with additional stuff attached as attributes.
makeObjectiveFunction = function(
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
  constraint.fn = NULL) {

  # sanity checks
  if (!is.null(name))
    assertString(name, na.ok = FALSE)
  if (!is.null(id))
    assertString(id, na.ok = FALSE)

  if (!is.null(description))
    assertString(description, na.ok = FALSE)
  assertFunction(fn)
  assertFlag(has.simple.signature, na.ok = FALSE)

  if (has.simple.signature) {
    fn = makeInternalObjectiveFunction(fn)
    if (!is.null(fn.mean)) {
      fn.mean = makeInternalObjectiveFunction(fn)
    }
  }

  assertClass(par.set, "ParamSet")
  assertInt(n.objectives, na.ok = FALSE, lower = 1L)
  assertFlag(noisy, na.ok = FALSE)
  if (!noisy && !is.null(fn.mean)) {
    stopf("Setting fn.mean only makes sense for noisy functions.")
  }
  if (!is.null(fn.mean)) {
    assertFunction(fn.mean)
  }
  assertLogical(minimize, len = n.objectives, any.missing = FALSE, all.missing = FALSE)
  assertFlag(vectorized, na.ok = FALSE)

  if (!has.simple.signature && vectorized) {
    stopf("At the moment we allow 'vectorized' functions only for functions with simple signature.")
  }

  if (!is.null(constraint.fn)) {
    assertFunction(constraint.fn)
  }

  structure(
    fn,
    name = coalesce(name, ""),
    id = coalesce(id, NA),
    description = coalesce(description, ""),
    par.set = par.set,
    noisy = noisy,
    fn.mean = fn.mean,
    minimize = minimize,
    vectorized = vectorized,
    constraint.fn = constraint.fn,
    n.objectives = n.objectives,
    class = c("smoof_function", "function")
  )
}
