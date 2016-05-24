#' Get the first/last element of a list/vector.
#'
#' @param x [\code{list} | \code{vector}]\cr
#'   The list or vector.
#' @return Selected element. The element name is dropped.
#' @export
getFirst = function(x) {
  assertVector(x, min.len = 1L)
  x[[1L]]
}

#' @rdname getFirst
#' @export
getLast = function(x) {
  assertVector(x, min.len = 1L)
  x[[length(x)]]
}
