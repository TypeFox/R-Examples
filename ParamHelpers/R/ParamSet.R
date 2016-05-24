#' @title Construct a parameter set.
#'
#' @description
#' \code{makeParamSet}: Contruct from a bunch of parameters.
#'
#' Multiple sets can be concatenated with \code{c}.
#'
#' The constructed S3 class is simply a list that contains the element \code{pars}.
#' \code{pars} is a list of the passed parameters, named by their ids.
#'
#' @param ... [\code{\link{Param}}]\cr
#'   Parameters.
#' @param params [list of \code{\link{Param}}]\cr
#'   List of parameters, alternative way instead of using \code{...}.
#' @param forbidden [\code{NULL} | R expression]\cr
#'   States forbidden region of parameter set via an expression.
#'   Every setting which satisfies this expression is considered to be infeasible.
#'   This makes it possible to exclude more complex region of the parameter space
#'   than through simple constraints or \code{requires}-conditions
#'   (although these should be always used when possible).
#'   If parameters have associated trafos, the forbidden region must always be specified on the original
#'   scale and not the transformed one.
#'   Default is \code{NULL} which means no forbidden region.
#' @return [\code{\link{ParamSet}}].
#' @aliases ParamSet
#' @export
#' @examples
#' makeParamSet(
#'   makeNumericParam("u", lower=1),
#'   makeIntegerParam("v", lower=1, upper=2),
#'   makeDiscreteParam("w", values=1:2),
#'   makeLogicalParam("x"),
#'   makeDiscreteVectorParam("y", len=2, values=c("a", "b"))
#' )
makeParamSet = function(..., params = NULL, forbidden = NULL) {
  pars = list(...)
  if (length(pars) > 0 && !is.null(params))
    stop("You can only use one of ... or params!")
  if (!is.null(params)) {
    assertList(params, types = "Param")
    pars = params
  } else {
    assertList(pars, types = "Param")
  }
  ns = extractSubList(pars, "id")
  if (any(duplicated(ns)))
    stop("All parameters must have unique names!")
  names(pars) = ns
  return(makeS3Obj("ParamSet", pars = pars, forbidden = forbidden))
}

getParSetPrintData = function(x, trafo = TRUE, used = TRUE, constr.clip = 40L) {
  d = lapply(x$pars, getParPrintData, trafo = trafo, used = used, constr.clip = constr.clip)
  return(do.call(rbind, d))
}

#' @export
print.ParamSet = function(x, ..., trafo = TRUE, used = TRUE, constr.clip = 40L) {
  if (isEmpty(x)) {
    print("Empty parameter set.")
  } else  {
    print(getParSetPrintData(x, trafo = trafo, used = used, constr.clip = constr.clip))
  }
  if (hasForbidden(x))
    catf("Forbidden region specified.")
  return(invisible(NULL))
}

#' @export
c.ParamSet = function(..., recursive = FALSE) {
  pss = list(...)
  pars = Reduce(c, lapply(pss, function(ps) ps$pars))
  # remove the names here. if 'params' is a par name, this wont work in the contructor call
  # but we are allowed to pass the list without names, as they are set again automatically later for pars
  names(pars) = NULL
  return(do.call(makeParamSet, pars))
}

#' Check whether parameter set is empty.
#'
#' @param par.set \code{\link{ParamSet}}]\cr
#'   Parameter set.
#' @return [\code{logical(1)}].
#' @export
isEmpty = function(par.set) {
  assertClass(par.set, "ParamSet")
  UseMethod("isEmpty")
}

#' @export
isEmpty.ParamSet = function(par.set) {
  return(length(par.set$pars) == 0)
}

#' \code{makeNumericParamSet}: Convenience function for numerics.
#'
#' @param id [\code{character(1)}]
#'   Name of parameter.
#' @param len [\code{integer(1)}]\cr
#'   Length of vector.
#' @param lower [\code{numeric}]\cr
#'   Lower bound.
#'   Default is \code{-Inf}.
#' @param upper [\code{numeric}] \cr
#'   Upper bound.
#'   Default is \code{Inf}.
#' @param vector [\code{logical(1)}] \cr
#'   Should a \code{NumericVectorParam} be used instead of
#'   n \code{NumericParam} objects?
#'   Default is \code{TRUE}.
#' @rdname makeParamSet
#' @export
makeNumericParamSet = function(id = "x", len, lower = -Inf, upper = Inf, vector = TRUE) {
  assertString(id)
  if (missing(len)) {
    if (!missing(lower))
      len = length(lower)
    else if (!missing(upper))
      len = length(upper)
  } else {
    len = asInt(len)
  }
  if (is.numeric(lower) && length(lower) == 1)
    lower = rep(lower, len)
  if (is.numeric(upper) && length(upper) == 1)
    upper = rep(upper, len)
    assertNumeric(lower, len = len)
    assertNumeric(upper, len = len)
    assertFlag(vector)
  if (vector) {
    return(makeParamSet(makeNumericVectorParam(id = id, len = len, lower = lower, upper = upper)))
  } else {
    return(makeParamSet(params = lapply(1:len, function(i)
      makeNumericParam(id = paste(id, i, sep = ""), lower = lower[i], upper = upper[i]))
    ))
  }
}


