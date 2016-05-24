#' Create a description object for a parameter.
#'
#' For each parameter type a special constructor function is available, see below.
#'
#' The S3 class is a list which stores these elements:
#' \describe{
#' \item{id [\code{character(1)}]}{See argument of same name.}
#' \item{type [\code{character(1)}]}{Data type of parameter. Possible types are \dQuote{numeric}, \dQuote{numericvector}, \dQuote{integer}, \dQuote{integervector}, \dQuote{logical}, \dQuote{logicalvector}, \dQuote{discrete}, \dQuote{discretevector}, \dQuote{function}, \dQuote{untyped}.}
#' \item{len [\code{integer(1)}]}{See argument of same name.}
#' \item{lower [\code{numeric}]}{See argument of same name. Length of this vector is \code{len}.}
#' \item{upper [\code{numeric}]}{See argument of same name. Length of this vector is \code{len}.}
#' \item{values [\code{list}]}{Discrete values, always stored as a named list.}
#' \item{trafo [\code{NULL} | \code{function(x)}]}{See argument of same name.}
#' \item{requires [\code{NULL} | \code{expression}]}{See argument of same name.}
#' }
#'
#' @param id [\code{character(1)}]\cr
#'   Name of parameter.
#' @param len [\code{integer(1)}]\cr
#'   Length of vector parameter.
#' @param lower [\code{numeric}]\cr
#'   Lower bounds.
#'   A singe value of length 1 is automatically replicated to \code{len} for vector parameters.
#'   Default is \code{-Inf}.
#' @param upper [\code{numeric}]\cr
#'   Upper bounds.
#'   A singe value of length 1 is automatically replicated to \code{len} for vector parameters.
#'   Default is \code{Inf}.
#' @param values [\code{vector} | \code{list}]\cr
#'   Possible discrete values. Instead of using a vector of atomic values,
#'   you are also allowed to pass a list of quite \dQuote{complex} R objects,
#'   which are used as discrete choices. If you do the latter,
#'   the elements must be uniquely named, so that the names can be used
#'   as internal representations for the choice.
#' @param cnames [\code{character}]\cr
#'   Component names for vector params (except discrete).
#'   Every function in this package that creates vector values for such a param, will name
#'   that vector with \code{cnames}.
#' @param allow.inf [\code{logical(1)}]\cr
#'   Allow infinite values for numeric and numericvector params to be feasible settings.
#'   Default is \code{FALSE}.
#' @param default [any]\cr
#'   Default value used in learner.
#'   If this argument is missing, it means no default value is available.
#' @param trafo [\code{NULL} | \code{function(x)}]\cr
#'   Function to transform parameter. It should be applied to the parameter value
#'   before it is, e.g., passed to a corresponding objective function.
#'   Function must accept a parameter value as the first argument and return a transformed one.
#'   Default is \code{NULL} which means no transformation.
#' @param requires [\code{NULL} | R expression]\cr
#'   States requirements on other paramaters' values, so that setting
#'   this parameter only makes sense if its requirements are satisfied (dependent parameter).
#'   Only really useful if the parameter is included in a \code{\link{ParamSet}}.
#'   Note that if your dependent parameter is a logical Boolean you need to verbosely write
#'   \code{requires = quote(a == TRUE)} and not \code{requires = quote(a)}.
#'   Default is \code{NULL} which means no requirements.
#' @param tunable [\code{logical(1)}]\cr
#'   Is this parameter tunable?
#'   Defining a parameter to be not-tunable allows to mark arguments like, e.g., \dQuote{verbose} or
#'   other purely technical stuff, and allows them to be excluded from later automatic optimization
#'   procedures that would try to consider all available parameters.
#'   Default is \code{TRUE} (except for \code{untyped}, \code{function}, \code{character} and
#'   \code{characterVector}) which means it is tunable.
#' @return [\code{\link{Param}}].
#' @name Param
#' @rdname Param
#' @examples
#' makeNumericParam("x",lower = -1, upper = 1)
#' makeNumericVectorParam("x", len = 2)
#' makeDiscreteParam("y", values = c("a","b"))
#' makeCharacterParam("z")
NULL

makeParam = function(id, type, len, lower, upper, values, cnames, allow.inf = FALSE, default,
  trafo = NULL, requires = NULL, tunable = TRUE) {

  #We cannot check default} for NULL or NA as this could be the default value!
  if (missing(default)) {
    has.default = FALSE
    default = NULL
  } else {
    has.default = TRUE
  }
  #FIXME: Do we need to check for NA here? hopefully not because this might occur in mlr?
  if (has.default && isScalarNA(default))
    warningf("NA used as a default value for learner parameter %s.\nParamHelpers uses NA as a special value for dependent parameters.", id)
  p = makeS3Obj("Param",
    id = id,
    type = type,
    len = len,
    lower = lower,
    upper = upper,
    values = values,
    cnames = cnames,
    allow.inf = allow.inf,
    has.default = has.default,
    default = default,
    trafo = trafo,
    requires = requires,
    tunable = tunable
  )
  if (has.default && !isFeasible(p, default))
    stop(p$id, " : 'default' must be a feasible parameter setting.")
  return(p)
}

getParPrintData = function(x, trafo = TRUE, used = TRUE, constr.clip = 40L) {
  g = function(n) collapse(sprintf("%.3g", n))
  if (isNumeric(x, include.int = TRUE))
    constr = sprintf("%s to %s", g(x$lower), g(x$upper))
  else if (isDiscrete(x, include.logical = FALSE))
    constr = clipString(collapse(names(x$values)), constr.clip)
  else
    constr = "-"
  d = data.frame(
    Type = x$type,
    len = ifelse(isVector(x), x$len, "-"),
    Def = if (x$has.default) paramValueToString(x, x$default) else "-",
    Constr = constr,
    Req = ifelse(is.null(x$requires), "-", "Y"),
    Tunable = x$tunable,
    stringsAsFactors = FALSE
  )
  if (trafo)
    d$Trafo = ifelse(is.null(x$trafo), "-", "Y")
  return(d)
}


#' @export
print.Param = function(x, ..., trafo = TRUE) {
  print(getParPrintData(x, trafo = trafo))
}

