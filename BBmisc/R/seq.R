#' Generate sequences along rows or cols.
#'
#' A simple convenience wrapper around \code{\link[base]{seq_len}}.
#'
#' @param x [\code{data.frame} | \code{matrix}]\cr
#'   Data frame, matrix or any object which supports \code{\link[base]{nrow}}
#'   or \code{\link[base]{ncol}}, respectively.
#' @return Vector of type [\code{integer}].
#' @export
#' @examples
#' data(iris)
#' seq_row(iris)
#' seq_col(iris)
seq_row = function(x) {
  seq_len(nrow(x))
}

#' @export seq_col
#' @rdname seq_row
seq_col = function(x) {
  seq_len(ncol(x))
}
