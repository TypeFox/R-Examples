#' @title Check if function should be minimized.
#'
#' @description Functions can have an associated global optimum. In this case
#' one needs to know whether the optimum is a minimum or a maximum.
#'
#' @template arg_smoof_function
#' @return [\code{logical}] Each component indicates whether the corresponding
#' objective should be minimized.
#' @export
shouldBeMinimized = function(fn) {
  UseMethod("shouldBeMinimized")
}

#' @export
shouldBeMinimized.smoof_function = function(fn) {
  return(attr(fn, "minimize"))
}

#' @export
shouldBeMinimized.smoof_wrapped_function = function(fn) {
  return(shouldBeMinimized(getWrappedFunction(fn)))
}
