#' Get parameter subset of only certain parameters.
#'
#' Parameter order is not changed.
#'
#' @template arg_parset
#' @param ids [\code{NULL} | \code{character}]\cr
#'   Vector with id strings containing the parameters to select. Has to be a
#'   subset of the parameter names within the parameter set.
#'   Per default (\code{ids = NULL}) no filtering based on names is done.
#' @param type [\code{NULL} | \code{character}]\cr
#'   Vector of allowed types, subset of: \dQuote{numeric}, \dQuote{integer}, \dQuote{numericvector},
#'   \dQuote{integervector}, \dQuote{discrete}, \dQuote{discretevector}, \dQuote{logical},
#'   \dQuote{logicalvector}, \dQuote{character}, \dQuote{charactervector},
#'   \dQuote{function}, \dQuote{untyped}.
#'   Setting \code{type = NULL}, which is the default, allows the consideration of all types.
#' @param tunable [\code{logical}]\cr
#'   Vector of allowed values for the property \code{tunable}. Accepted arguments are
#'   \code{TRUE}, \code{FALSE} or \code{c(TRUE, FALSE)}.
#'   The default is \code{c(TRUE, FALSE)}, i.e. none of the parameters will be filtered out.
#' @return [\code{\link{ParamSet}}].
#' @examples
#' ps = makeParamSet(
#'   makeNumericParam("u", lower = 1),
#'   makeIntegerParam("v", lower = 1, upper = 2),
#'   makeDiscreteParam("w", values = 1:2),
#'   makeLogicalParam("x"),
#'   makeCharacterParam("s"),
#'   makeNumericParam("y", tunable = FALSE)
#' )
#'
#' # filter for numeric and integer parameters
#' filterParams(ps, type = c("integer", "numeric"))
#'
#' # filter for tunable, numeric parameters
#' filterParams(ps, type = "numeric", tunable = TRUE)
#'
#' # filter for all numeric parameters among "u", "v" and "x"
#' filterParams(ps, type = "numeric", ids = c("u", "v", "x"))
#' @export
filterParams = function(par.set, ids = NULL, type = NULL, tunable = c(TRUE, FALSE)) {
  # FIXME: how do we handle this, this also affects "requires" the same way?
  # if we drop same params the expressions can become invalid?
  # if (!is.null(par.set$forbidden))
    # stopf("Operation not allowed for param set with forbidden region currently!")
  if (!is.null(ids)) {
    assertSubset(ids, names(par.set$pars))
    par.set$pars = Filter(function(p) p$id %in% ids, par.set$pars)
  }
  if (!is.null(type)) {
    assertSubset(type, c("numeric", "integer", "numericvector", "integervector", "discrete",
        "discretevector", "logical", "logicalvector", "character", "charactervector",
        "function", "untyped"))
    par.set$pars = Filter(function(p) p$type %in% type, par.set$pars)
  }
  assertLogical(tunable, min.len = 1L, max.len = 2L, unique = TRUE)
  par.set$pars = Filter(function(p) p$tunable %in% tunable, par.set$pars)
  return(par.set)
}
