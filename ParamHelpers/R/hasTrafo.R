#' Check parameter / parameter set for trafos.
#'
#' \code{TRUE} iff the parameter has any trafos or any parameter in the set has
#' trafos.
#'
#' @template arg_par_or_set
#' @return [\code{logical(1)}].
#' @export
hasTrafo = function(par) {
  UseMethod("hasTrafo")
}

#' @export
hasTrafo.Param = function(par) {
  return(!is.null(par$trafo))
}

#' @export
hasTrafo.ParamSet = function(par) {
  return(any(vapply(par$pars, hasTrafo, logical(1L))))
}

