#' Return the ID / short name of the function.
#'
#' @template arg_smoof_function
#' @return [\code{character(1)}]
#' @export
getID = function(fn) {
  UseMethod("getID")
}

#' @export
getID.smoof_function = function(fn) {
  return(attr(fn, "id"))
}

#' @export
getID.smoof_wrapped_function = function(fn) {
  return(getID(getWrappedFunction(fn)))
}
