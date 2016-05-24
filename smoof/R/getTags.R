#' Returns vector of associated tags.
#'
#' @template arg_smoof_function
#' @return [\code{character}]
#' @export
getTags = function(fn) {
  UseMethod("getTags")
}

#' @export
getTags.smoof_function = function(fn) {
  return(attr(fn, "tags"))
}

#' @export
getTags.smoof_wrapped_function = function(fn) {
  return(getTags(getWrappedFunction(fn)))
}

#' @export
getTags.smoof_generator = function(fn) {
  return(attr(fn, "tags"))
}
