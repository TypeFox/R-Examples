#' @title Check if parameter value is valid.
#'
#' @description
#' Check if a parameter value satisfies the constraints of the
#' parameter description. This includes the \code{requires} expressions and
#' the \code{forbidden} expression, if \code{par} is a \code{\link{ParamSet}}.
#' If \code{requires} is not satisfied,
#' the parameter value must be set to scalar \code{NA} to be still feasible, a single scalar even in a
#' case of a vector parameter.
#'
#' If the parameter has \code{cnames}, these are also checked.
#'
#' @template arg_par_or_set
#' @param x [any] \cr
#'   Single value to check.
#'   For a parameter set this must be a list.
#'   If the list is named, it is possible to only pass a subset of parameters defined
#'   in the \code{\link{ParamSet}} \code{par}. In that case, only conditions regarding the passed
#'   parameters are checked.
#'   (Note that this might not work if one of the passed params has a \code{requires} setting
#'   which refers to an unpassed param.)
#' @return logical(1)
#' @examples
#' p = makeNumericParam("x", lower = -1, upper = 1)
#' isFeasible(p, 0) # True
#' isFeasible(p, 2) # False, out of bounds
#' isFeasible(p, "a") # False, wrong type
#' # now for parameter sets
#' ps = makeParamSet(
#'   makeNumericParam("x", lower = -1, upper = 1),
#'   makeDiscreteParam("y", values = c("a", "b"))
#' )
#' isFeasible(ps, list(0, "a")) # True
#' isFeasible(ps, list("a", 0)) # False, wrong order
#' @export
isFeasible = function(par, x) {
  UseMethod("isFeasible")
}

#' @export
isFeasible.Param = function(par, x) {
  # we dont have to consider requires here, it is not a param set
  constraintsOkParam(par, x)
}

#' @export
isFeasible.LearnerParam = function(par, x) {
  # we dont have to consider requires here, it is not a param set
  constraintsOkLearnerParam(par, x)
}

#' @export
isFeasible.ParamSet = function(par, x) {
  named = testNamed(x)
  if (!is.list(x) || !(!named || all(names(x) %in% names(par$pars)) || !(named || length(x) == length(par$pars))))
    return(FALSE)
  if (isForbidden(par, x))
    return(FALSE)
  if (named) {
    par = filterParams(par, ids = names(x))
    x = x[names(par$pars)]
  }
  #FIXME: very slow
  for (i in seq_along(par$pars)) {
    p = par$pars[[i]]
    v = x[[i]]
    # no requires, just check constraints
    if(is.null(p$requires)) {
      if (!isFeasible(p, v))
        return(FALSE)
    } else {
      # requires, is it ok?
      if (!requiresOk(par, x, i)) {
        # if not, val must be NA
        if (!isScalarNA(v))
          return(FALSE)
      } else {
        # requires, ok, check constraints
        if (!isFeasible(p, v))
          return(FALSE)
      }
    }
  }
  return(TRUE)
}

# are the contraints ok for value of a param (not considering requires)
constraintsOkParam = function(par, x) {
  type = par$type
  # this should work in any! case.
  if (type == "untyped")
    return(TRUE)
  inValues = function(v) any(sapply(par$values, function(w) isTRUE(all.equal(w, v))))
  ok = if (type == "numeric")
    is.numeric(x) && length(x) == 1 && (par$allow.inf || is.finite(x)) && x >= par$lower && x <= par$upper
  else if (type == "integer")
    is.numeric(x) && length(x) == 1 && is.finite(x) && x >= par$lower && x <= par$upper && x == as.integer(x)
  else if (type == "numericvector")
    is.numeric(x) && length(x) == par$len && all((par$allow.inf | is.finite(x)) & x >= par$lower & x <= par$upper)
  else if (type == "integervector")
    is.numeric(x) && length(x) == par$len && all(is.finite(x) & x >= par$lower & x <= par$upper & x == as.integer(x))
  else if (type == "discrete")
    inValues(x)
  else if (type == "discretevector")
    is.list(x) && length(x) == par$len && all(sapply(x, inValues))
  else if (type == "logical")
    is.logical(x) && length(x) == 1 && !is.na(x)
  else if (type == "logicalvector")
    is.logical(x) && length(x) == par$len && !any(is.na(x))
  else if (type == "character")
    is.character(x) && length(x) == 1 && !is.na(x)
  else if (type == "charactervector")
    is.character(x) && length(x) == par$len && !any(is.na(x))
  else if (type == "function")
    is.function(x)

  # if we have cnames, check them
  if (!is.null(par$cnames))
    ok = ok && !is.null(names(x)) && all(names(x) == par$cnames)

  return(ok)
}

constraintsOkLearnerParam = function(par, x) {
  inValues = function(v) any(sapply(par$values, function(w) isTRUE(all.equal(w, v))))
  type = par$type
  # extra case for unkown dim in vector
  if (type == "numericvector")
    is.numeric(x) && (is.na(par$len) || length(x) == par$len) && all(is.finite(x) & x >= par$lower & x <= par$upper)
  else if (type == "integervector")
    is.numeric(x) && (is.na(par$len) || length(x) == par$len) && all(is.finite(x) & x >= par$lower & x <= par$upper & x == as.integer(x))
  else if (type == "logicalvector")
    is.logical(x) && (is.na(par$len) || length(x) == par$len) && !any(is.na(x))
  else if (type == "discretevector")
    is.list(x) && (is.na(par$len) || length(x) == par$len) && all(sapply(x, inValues))
  else
    constraintsOkParam(par, x)
}

# is the requires part of the ith param valid for value x (x[[i]] is value or ith param)
# assumes that param actually has a requires part
requiresOk = function(par.set, x, i) {
  if (is.null(par.set$pars[[i]]$requires)) {
    TRUE
  } else {
    isTRUE(eval(par.set$pars[[i]]$requires, envir = x))
  }
}
