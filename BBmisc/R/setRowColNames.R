#' Wrapper for \code{rownames(x) = y}, \code{colnames(x) = y}.
#'
#' @param x [\code{matrix} | \code{data.frame}]\cr
#'   Matrix or data.frame.
#' @param names [\code{character}]\cr
#'   New names for rows / columns.  
#' @return Changed object \code{x}.
#' @export
#' @examples
#' setColNames(matrix(1:4, 2, 2), c("a", "b"))
setRowNames = function(x, names) {
  rownames(x) = names
  x
}

#' @rdname setRowNames
#' @export
setColNames = function(x, names) {
  colnames(x) = names
  x
}
