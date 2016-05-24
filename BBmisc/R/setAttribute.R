#' A wrapper for \code{attr(x, which) = y}.
#'
#' @param x [any]\cr
#'  Your object.
#' @param which [\code{character(1)}]\cr
#'  Name of the attribute to set
#' @param value [\code{ANY}]\cr
#'  Value for the attribute.
#' @return Changed object \code{x}.
#' @export
#' @examples
#' setAttribute(list(), "foo", 1)
setAttribute = function(x, which, value) {
  attr(x, which) = value
  x
}
