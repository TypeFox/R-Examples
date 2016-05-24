#' A wrapper for \code{class(x) = classes}.
#'
#' @param x [any]\cr
#'   Your object.
#' @param classes [\code{character}]\cr
#'  New classes.
#' @return Changed object \code{x}.
#' @export
#' @examples
#' setClasses(list(), c("foo1", "foo2"))
setClasses = function(x, classes) {
  class(x) = classes
  x
}
