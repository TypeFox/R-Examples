#' Extracts a named element from a list of lists.
#'
#' @param xs [\code{list}]\cr
#'   A list of vectors of the same length.
#' @param row.names [\code{character} | \code{integer} | \code{NULL}]\cr
#'   Row names of result.
#'   Default is to take the names of the elements of \code{xs}.
#' @param col.names [\code{character} | \code{integer} | \code{NULL}]\cr
#'   Column names of result.
#'   Default is to take the names of the elements of \code{xs}.
#' @return [\code{matrix}].
#' @export
asMatrixCols = function(xs, row.names, col.names) {
  assertList(xs)
  n = length(xs)
  if (n == 0L)
    return(matrix(0, nrow = 0L, ncol = 0L))
  assertList(xs, types = "vector")

  m = unique(viapply(xs, length))
  if (length(m) != 1L)
    stopf("Vectors must all be of the same length!")

  if (missing(row.names)) {
    row.names = names(xs[[1L]])
  }

  if (missing(col.names)) {
    col.names = names(xs)
  }

  xs = unlist(xs)
  dim(xs) = c(m, n)
  rownames(xs) = row.names
  colnames(xs) = col.names
  return(xs)
}

#' @rdname asMatrixCols
#' @export
asMatrixRows = function(xs, row.names, col.names) {
  t(asMatrixCols(xs, row.names = col.names, col.names = row.names))
}
