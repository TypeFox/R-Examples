#' Return the name of the function.
#'
#' @template arg_smoof_function
#' @return [\code{character(1)}]
#' @export
getName = function(fn) {
  UseMethod("getName")
}

#' @export
getName.smoof_function = function(fn) {
  return(attr(fn, "name"))
}

#' @export
getName.smoof_wrapped_function = function(fn) {
  return(getName(getWrappedFunction(fn)))
}
