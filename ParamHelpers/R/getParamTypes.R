#' Returns type information for a parameter set.
#'
#' @template arg_parset
#' @param df.cols [\code{logical(1)}]\cr
#'   If \code{FALSE} simply return the parameter types in the set,
#'   i.e., \code{par$type}.
#'   If \code{TRUE}, convert types so they correspond to R types of a data.frame
#'   where values of this set might be stored.
#'   This also results in replication of output types for
#'   vector parameters.
#'   Default is \code{FALSE}.
#' @param df.discretes.as.factor [\code{logical(1)}]\cr
#'   If \code{df.cols} is \code{TRUE}:
#'   Should type for discrete params be \code{factor} or \code{character}?
#'   Default is \code{TRUE}.
#' @param use.names [\code{logical(1)}]\cr
#'   Name the result vector?
#'   Default is \code{FALSE}.
#' @param with.nr [\code{logical(1)}]\cr
#'   Should number from 1 to length be appended to name?
#'   Only used if \code{use.names} and \code{df.cols} are \code{TRUE}.
#'   Default is \code{TRUE}.
#' @return [\code{character}].
#' @export
getParamTypes = function(par.set, df.cols = FALSE, df.discretes.as.factor = TRUE,
  use.names = FALSE, with.nr = TRUE) {
  assertClass(par.set, "ParamSet")
  assertFlag(df.cols)
  assertFlag(df.discretes.as.factor)
  assertFlag(use.names)
  assertFlag(with.nr)

  types = extractSubList(par.set$pars, "type")
  recode = function(types, ...) {
    args = as.character(list(...))
    for (i in seq(1, length(args), 2)) {
      types[types == args[i]] = args[i + 1]
    }
    types = rep(types, getParamLengths(par.set))
    return(types)
  }

  if (df.cols) {
    types = if (df.discretes.as.factor) {
      recode(types,
        "numericvector", "numeric",
        "integervector", "integer",
        "discrete", "factor",
        "discretevector", "factor",
        "logicalvector", "logical",
        "charactervector", "character"
      )
    } else {
      recode(types,
        "numericvector", "numeric",
        "integervector", "integer",
        "discrete", "character",
        "discretevector", "character",
        "logicalvector", "logical",
        "charactervector", "character"
      )
    }
  }

  ns = if (use.names)
    getParamIds(par.set, repeated = df.cols, with.nr = with.nr)
  else
    NULL
  return(setNames(types, ns))
}
