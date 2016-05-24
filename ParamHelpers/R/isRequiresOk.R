#' @title Check if parameter requirements are met.
#'
#' @description
#' Check if a parameter value satisfies the requirements of the
#' parameter description. This only checks the \code{requires} expressions.
#'
#' @template arg_parset
#' @param par.vals \code{list()} \cr
#'   List of parameter settings.
#' @param ids \code{character()} \cr
#'   \code{id}s of the param.vals to check. Default is \code{names(par.vals)}.
#' @param use.defaults \code{logical()} \cr
#'   Some requirements relay on default values of the \code{par.set}. Default is \code{TRUE}, which means that if the value is not present in \code{par.vals} the default value will be considered.
#' @return \code{logical(1)} \cr
#' @export
isRequiresOk = function(par.set, par.vals, ids = names(par.vals), use.defaults = TRUE) {
  assertClass(par.set, "ParamSet")
  assertList(par.vals)
  assertNamed(par.vals)
  if (is.numeric(ids))
    assertInteger(ids, lower = 1L, upper = length(par.vals), unique = TRUE)
  else
    assertSubset(ids, choices = names(par.vals))
  assertFlag(use.defaults)
  if (use.defaults) {
    par.vals.env = insert(getDefaults(par.set), par.vals)
  } else {
    par.vals.env = par.vals
  }
  requireOks = vlapply(names(par.vals), function(par.name) {
    requiresOk(par.set, par.vals.env, par.name)
  })
  if (any(!requireOks)) {
    #just constructing an informative error message
    #FIXME Maybe use paramValueToString
    par.names.failed = names(requireOks)[!requireOks]
    par.vals.failed = par.vals[par.names.failed]
    requires.failed = as.character(extractSubList(par.set$pars, "requires")[par.names.failed])
    stopf("The following parameter settings do not meet the requirements: %s",
      paste0(par.names.failed, "=", par.vals.failed, " needs ", requires.failed, collapse = ", "))
  }
  return(all(requireOks))
}
