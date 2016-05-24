#' Are all elements of a list / vector uniquely named?
#'
#' \code{NA} or \dQuote{} are not allowed as names.
#'
#' @param x [\code{vector}]\cr
#'   The vector or list.
#' @return [\code{logical(1)}].
#' @export
#' @examples
#' isProperlyNamed(list(1))
#' isProperlyNamed(list(a = 1))
#' isProperlyNamed(list(a = 1, 2))
isProperlyNamed = function(x) {
  ns = names2(x)
  length(x) == 0L || !(any(is.na(ns)) || anyDuplicated(ns))
}
