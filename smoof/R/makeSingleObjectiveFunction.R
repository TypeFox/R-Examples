#' Generator for single-objective target functions.
#'
#' @template arg_name
#' @template arg_id
#' @template arg_description
#' @template arg_fn
#' @template arg_has_simple_signature
#' @template arg_par_set
#' @template arg_noisy
#' @template arg_fn_mean
#' @template arg_minimize
#' @template arg_vectorized
#' @template arg_constraint_fn
#' @param tags [\code{character}]\cr
#'   Optional character vector of tags or keywords which characterize the function,
#'   e.~g. \dQuote{unimodal}, \dQuote{separable}. See \code{\link{getAvailableTags}} for
#'   a character vector of allowed tags.
#' @param global.opt.params [\code{list} | \code{numeric} | \code{data.frame} | \code{matrix} | \code{NULL}]\cr
#'   Default is \code{NULL} which means unknown. Passing a \code{numeric} vector will
#'   be the most frequent case (numeric only functions). In this case there is only a
#'   single global optimum. If there are multiple global optima, passing a numeric
#'   \code{matrix} is the best choice. Passing a \code{list} or a \code{data.frame}
#'   is necessary if your function is mixed, e.g., it expects both numeric and discrete
#'   parameters. Internally, however, each representation is casted to a \code{data.frame}
#'   for reasons of consistency.
#' @param global.opt.value [\code{numeric(1)} | \code{NULL}]\cr
#'   Global optimum value if known. Default is \code{NULL}, which means unknown. If
#'   only the \code{global.opt.params} are passed, the value is computed automatically.
#' @return [\code{function}] Objective function with additional stuff attached as attributes.
#' @examples
#' library(ggplot2)
#'
#' fn = makeSingleObjectiveFunction(
#'   name = "Sphere Function",
#'   fn = function(x) sum(x^2),
#'   par.set = makeNumericParamSet("x", len = 1L, lower = -5L, upper = 5L),
#'   global.opt.params = list(x = 0)
#' )
#' print(fn)
#' print(autoplot(fn))
#'
#' fn.num2 = makeSingleObjectiveFunction(
#'   name = "Numeric 2D",
#'   fn = function(x) sum(x^2),
#'   par.set = makeParamSet(
#'     makeNumericParam("x1", lower = -5, upper = 5),
#'     makeNumericParam("x2", lower = -10, upper = 20)
#'   )
#' )
#' print(fn)
#' print(autoplot(fn))
#'
#' fn.mixed = makeSingleObjectiveFunction(
#'   name = "Mixed 2D",
#'   fn = function(x) x$num1^2 + as.integer(as.character(x$disc1) == "a"),
#'   has.simple.signature = FALSE,
#'   par.set = makeParamSet(
#'     makeNumericParam("num1", lower = -5, upper = 5),
#'     makeDiscreteParam("disc1", values = c("a", "b"))
#'   ),
#'   global.opt.params = list(num1 = 0, disc1 = "b")
#' )
#' print(fn)
#' print(autoplot(fn))
#' @export
makeSingleObjectiveFunction = function(
  name = NULL,
  id = NULL,
  description = NULL,
  fn,
  has.simple.signature = TRUE,
  vectorized = FALSE,
  par.set,
  noisy = FALSE,
  fn.mean = NULL,
  minimize = TRUE,
  constraint.fn = NULL,
  tags = character(0),
  global.opt.params = NULL,
  global.opt.value = NULL) {

  smoof.fn = makeObjectiveFunction(
    name, id, description, fn,
    has.simple.signature, par.set, 1L,
    noisy, fn.mean, minimize, vectorized, constraint.fn
  )
  n.params = getNumberOfParameters(smoof.fn)

  #FIXME: currently we offer this only for single objective functions
  assertSubset(tags, choices = getAvailableTags(), empty.ok = TRUE)

  if (!is.null(global.opt.params)) {
    if (!testDataFrame(global.opt.params)) {
      # single numeric only value passed
      if (testNumeric(global.opt.params, len = n.params, any.missing = FALSE)) {
        global.opt.params = as.data.frame(t(global.opt.params))
      } else if (testList(global.opt.params, len = n.params, any.missing = FALSE)) {
        global.opt.params = as.data.frame(global.opt.params)
      } else if (testMatrix(global.opt.params)) {
        global.opt.params = as.data.frame(global.opt.params)
      } else {
        stopf("Parameter(s) for known global optima must be passed as vector, list, matrix or data.frame.")
      }
      colnames(global.opt.params) = getParamIds(par.set, with.nr = TRUE, repeated = TRUE)
    }
    assertDataFrame(global.opt.params, ncols = n.params, col.names = "unique")

    # check if the passed parameters are indeed within the feasible region
    lapply(1:nrow(global.opt.params), function(i) {
      if (!isFeasible(par.set, ParamHelpers::dfRowToList(global.opt.params, par.set, i))) {
        stopf("Global optimum out of bounds.")
      }
    })
    if (!setequal(getParamIds(par.set, repeated = TRUE, with.nr = TRUE), colnames(global.opt.params))) {
      stopf("Names of values and parameter names do not match.")
    }
  }
  if (is.null(global.opt.value) && !is.null(global.opt.params)) {
    global.opt.value = smoof.fn(global.opt.params[1, ])
    assertNumber(global.opt.value, na.ok = FALSE, finite = TRUE)
  }

  smoof.fn = setAttribute(smoof.fn, "global.opt.params", global.opt.params)
  smoof.fn = setAttribute(smoof.fn, "global.opt.value", global.opt.value)
  smoof.fn = setAttribute(smoof.fn, "tags", tags)

  class(smoof.fn) = c("smoof_single_objective_function", class(smoof.fn))

  return(smoof.fn)
}

#' @export
print.smoof_function = function(x, ...) {
  n.objectives.text = ifelse(isSingleobjective(x), "Single", "Multi")
  catf("%s-objective function", n.objectives.text)
  if (isMultiobjective(x)) {
    catf("Number of objectives: %i", getNumberOfObjectives(x))
    ref.point = getRefPoint(x)
    if (!is.null(x)) {
      catf("Reference point:      (%s)", collapse(ref.point, ", "))
    }
  }
  catf("Name: %s", getName(x))
  description = getDescription(x)
  catf("Description: %s", if (description == "") "no description" else description)

  catf("Tags: %s", collapse(getTags(x), sep = ", "))
  catf("Noisy: %s", as.character(isNoisy(x)))
  catf("Minimize: %s", collapse(shouldBeMinimized(x)))
  catf("Constraints: %s", as.character(hasConstraints(x)))
  catf("Number of parameters: %i", getNumberOfParameters(x))
  print(getParamSet(x))
  if (hasGlobalOptimum(x)) {
    opt = getGlobalOptimum(x)
    catf("Global optimum objective value of %.4f at", opt$value)
    print(opt$param)
  }
}
