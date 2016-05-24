
#' A wrapper to add to the class attribute.
#'
#' @param x [any]\cr
#'   Your object.
#' @param classes [\code{character}]\cr
#'  Classes to add. Will be added in front (specialization).
#' @return Changed object \code{x}.
#' @export
#' @examples
#' addClasses(list(), c("foo1", "foo2"))
addClasses = function(x, classes) {
  class(x) = c(classes, class(x))
  x
}
