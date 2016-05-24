#' Return the description of the function.
#'
#' @template arg_smoof_function
#' @return [\code{character(1)}]
#' @export
getDescription = function(fn) {
  UseMethod("getDescription")
}

#' @export
getDescription.smoof_function = function(fn) {
  return(attr(fn, "description"))
}

#' @export
getDescription.smoof_wrapped_function = function(fn) {
  return(getDescription(getWrappedFunction(fn)))
}
