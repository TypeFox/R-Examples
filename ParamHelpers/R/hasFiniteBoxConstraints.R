#' Checks if a parameter or each parameter of a parameter set has ONLY
#' finite lower and upper bounds.
#'
#' @template arg_par_or_set
#' @return [\code{logical(1)}]
#' @export
hasFiniteBoxConstraints = function(par) {
  UseMethod("hasFiniteBoxConstraints")
}

#' @export
hasFiniteBoxConstraints.Param = function(par) {
  bounds = c(par$lower, par$upper)
  if (length(bounds) == 0)
    return(TRUE)
  return(all(is.finite(bounds)))
}

#' @export
hasFiniteBoxConstraints.ParamSet = function(par) {
  bounds = c(getLower(par), getUpper(par))
  if (length(bounds) == 0)
    return(TRUE)
  return(all(is.finite(bounds)))
}
