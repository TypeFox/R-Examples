#' Check parameter / parameter set for requirements / dependencies.
#'
#' \code{TRUE} iff the parameter has any requirements or any parameter in the set has
#' requirements.
#'
#' @template arg_par_or_set
#' @return [\code{logical(1)}].
#' @export
hasRequires = function(par) {
  UseMethod("hasRequires")
}

#' @export
hasRequires.Param = function(par) {
  return(!is.null(par$requires))
}

#' @export
hasRequires.ParamSet = function(par) {
  return(any(vapply(par$pars, hasRequires, logical(1L))))
}

