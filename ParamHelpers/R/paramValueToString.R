#' Convert a value to a string.
#'
#' Useful helper for logging.
#' For discrete parameter values always the name of the discrete value is used.
#'
#' @template arg_par_or_set
#' @param x [any]\cr
#'   Value for parameter or value for parameter set. In the latter case it must be named list.
#'   For discrete parameters their values must be used, not their names.
#' @param show.missing.values [\code{logical(1)}]\cr
#'   Display \dQuote{NA} for parameters, which have no setting, because their requirements are not
#'   satified (dependent parameters), instead of displaying nothing?
#'   Default is \code{FALSE}.
#' @param num.format [\code{character(1)}]\cr
#'   Number format for output of numeric parameters. See the details section of the manual for
#'   \code{\link[base]{sprintf}} for details.
#' @return [\code{character(1)}].
#' @export
#' @examples
#' p = makeNumericParam("x")
#' paramValueToString(p, 1)
#' paramValueToString(p, 1.2345)
#' paramValueToString(p, 0.000039)
#' paramValueToString(p, 8.13402, num.format = "%.2f")
#'
#' p = makeIntegerVectorParam("x", len=2)
#' paramValueToString(p, c(1L, 2L))
#'
#' p = makeLogicalParam("x")
#' paramValueToString(p, TRUE)
#'
#' p = makeDiscreteParam("x", values=list(a=NULL, b=2))
#' paramValueToString(p, NULL)
#'
#' ps = makeParamSet(
#'   makeNumericVectorParam("x", len=2L),
#'   makeDiscreteParam("y", values=list(a=NULL, b=2))
#' )
#' paramValueToString(ps, list(x=c(1,2), y=NULL))
paramValueToString = function(par, x, show.missing.values = FALSE, num.format = "%.3g") {
  assertFlag(show.missing.values)
  assertString(num.format)
  UseMethod("paramValueToString")
}

#' @export
paramValueToString.Param = function(par, x, show.missing.values = FALSE, num.format = "%.3g") {
  # handle missings
  if (isMissingValue(x)) {
    if (show.missing.values)
      return("NA")
    else
      return("")
  }

  type = par$type
  if (type == "numeric")
    sprintf(num.format, x)
  else if (type == "numericvector")
    paste(sprintf(num.format, x), collapse=",")
  else if (type == "integer")
    as.character(x)
  else if (type == "integervector")
    paste(as.character(x), collapse=",")
  else if (type == "logical")
    as.character(x)
  else if (type == "logicalvector")
    paste(as.character(x), collapse=",")
  else if (type == "discrete")
    discreteValueToName(par, x)
  else if (type == "discretevector")
    collapse(discreteValueToName(par, x))
  else if (type == "character")
    as.character(x)
  else if (type == "charactervector")
    collapse(x)
  else if (type == "function")
    "<function>"
  else if (type == "untyped")
    sprintf("<%s>", class(x)[1])
}

#' @export
paramValueToString.ParamSet = function(par, x, show.missing.values = FALSE, num.format = "%.3g") {
  assertList(x)
  if (!isProperlyNamed(x))
    stop("'x' must be a properly named list!")
  rest = setdiff(names(x), names(par$pars))
  if (length(rest) > 0L)
    stopf("Not all names of 'x' occur in par set 'par': %s", collapse(rest))
  res = character(0)
  for (i in seq_along(x)) {
    pn = names(x)[i]
    val = x[[pn]]
    if (show.missing.values || !isMissingValue(val))  {
      p = par$pars[[pn]]
      res[length(res)+1] = sprintf("%s=%s", pn, paramValueToString(p, val, show.missing.values, num.format))
    }
  }
  return(collapse(res, sep = "; "))
}
