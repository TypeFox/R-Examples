#' Convert rows (columns) of data.frame or matrix to lists.
#'
#' For each row, one list/vector is constructed, each entry of
#' the row becomes a list/vector element.
#'
#' @param x [\code{matrix} | \code{data.frame}]\cr
#'   Object to convert.
#' @param name.list [\code{logical(1)}]\cr
#'   Name resulting list with names of rows (cols) of \code{x}?
#'   Default is \code{FALSE}.
#' @param name.vector [\code{logical(1)}]\cr
#'   Name vector elements in resulting list with names of cols (rows) of \code{x}?
#'   Default is \code{FALSE}.
#' @param factors.as.char [\code{logical(1)}]\cr
#'   If \code{x} is a data.frame, convert factor columns to
#'   string elements in the resulting lists?
#'   Default is \code{TRUE}.
#' @param as.vector [\code{logical(1)}]\cr
#'   If \code{x} is a matrix, store rows as vectors in the resulting list - or otherwise as lists?
#'   Default is \code{TRUE}.
#' @return [\code{list} of lists or vectors].
#' @export
convertRowsToList = function(x, name.list = TRUE, name.vector = FALSE,
  factors.as.char = TRUE, as.vector = TRUE) {
  assert(checkMatrix(x), checkDataFrame(x))
  assertFlag(name.list)
  assertFlag(name.vector)
  assertFlag(factors.as.char)
  assertFlag(as.vector)
  ns.list = if (name.list) rownames(x) else NULL
  ns.vector = if (name.vector) colnames(x) else NULL
  if (is.matrix(x)) {
    if (as.vector)
      res = lapply(seq_row(x), function(i) setNames(x[i, ], ns.vector))
    else
      res = lapply(seq_row(x), function(i) setNames(as.list(x[i, ]), ns.vector))
  } else if (is.data.frame(x)) {
    if (factors.as.char)
      x = convertDataFrameCols(x, factors.as.char = TRUE)
    res = rowLapply(x, function(row) setNames(as.list(row), ns.vector))
  }
  setNames(res, ns.list)
}

#' @rdname convertRowsToList
#' @export
convertColsToList = function(x, name.list = FALSE, name.vector= FALSE,
  factors.as.char = TRUE, as.vector = TRUE) {

  # we need a special case for df and can ignore as.vector in it
  if (is.data.frame(x)) {
    if (factors.as.char)
      x = convertDataFrameCols(x, factors.as.char = TRUE)
    y = as.list(x)
    if (name.vector) {
      ns.vector = if (name.vector) colnames(x) else NULL
      y = lapply(y, function(z) setNames(z, ns.vector))
    }
    colnames(y) = if (name.list) colnames(x) else NULL
    return(y)
  }

  convertRowsToList(t(x), name.list = name.list, name.vector = name.vector,
    factors.as.char = factors.as.char, as.vector = as.vector)
}
