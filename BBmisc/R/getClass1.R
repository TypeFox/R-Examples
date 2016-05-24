#' Wrapper for \code{class(x)[1]}.
#'
#' @param x [any]\cr
#'   Input object.
#' @return [\code{character(1)}].
#' @note \code{getClass} is a function in \code{methods}. Do not confuse.
#' @export
getClass1 = function(x) {
  class(x)[1L]
}
