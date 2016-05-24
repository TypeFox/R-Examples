#' Apply function to rows of a data frame.
#'
#' Just like an \code{\link[base]{lapply}} on data frames,
#' but on the rows.
#'
#' @param df [\code{data.frame}]\cr
#'   Data frame.
#' @param fun [\code{function}]\cr
#'   Function to apply. Rows are passed as list or vector,
#'   depending on argument \code{unlist}, as first argument.
#' @param ... [\code{ANY}]\cr
#'   Additional arguments for \code{fun}.
#' @param unlist [\code{logical(1)}]\cr
#'   Unlist the row? Note that automatic conversion may be triggered for
#'   lists of mixed data types
#'   Default is \code{FALSE}.
#' @param simplify [\code{logical(1)} | character(1)]\cr
#'   Should the result be simplified?
#'   See \code{\link{sapply}}.
#'   If \dQuote{cols}, we expect the call results to be vectors of the same length and they are
#'   arranged as the columns of the resulting matrix.
#'   If \dQuote{rows}, likewise, but rows of the resulting matrix.
#'   Default is \code{TRUE}.
#' @param use.names [\code{logical(1)}]\cr
#'   Should result be named by the row names of \code{df}?
#'   Default is \code{TRUE}.
#' @return [\code{list} or simplified object]. Length is \code{nrow(df)}.
#' @export
#' @examples
#'  rowLapply(iris, function(x) x$Sepal.Length + x$Sepal.Width)
rowLapply = function(df, fun, ..., unlist = FALSE) {
  assertDataFrame(df)
  fun = match.fun(fun)
  assertFlag(unlist)
  if (unlist) {
    .wrap = function(.i, .df, .fun, ...)
      .fun(unlist(.df[.i, , drop = FALSE], recursive = FALSE, use.names = TRUE), ...)
  } else {
    .wrap = function(.i, .df, .fun, ...)
      .fun(as.list(.df[.i, , drop = FALSE]), ...)
  }

  lapply(seq_row(df), .wrap, .fun = fun, .df = df, ...)
}

#' @export
#' @rdname rowLapply
rowSapply = function(df, fun, ..., unlist = FALSE, simplify = TRUE, use.names = TRUE) {
  assert(checkFlag(simplify), checkChoice(simplify, c("cols", "rows")))
  assertFlag(use.names)
  ys = rowLapply(df, fun, ..., unlist = unlist)

  # simplify result
  if (length(ys) > 0L) {
    if (isTRUE(simplify)) {
      ys = simplify2array(ys)
    } else if (simplify == "rows") {
      ys = asMatrixRows(ys)
    } else if (simplify == "cols") {
      ys = asMatrixCols(ys)
    }
  }

  # set names
  if (use.names) {
    if (is.matrix(ys)) {
      colnames(ys) = rownames(df)
      rownames(ys) = NULL
    } else {
      names(ys) = rownames(df)
    }
  } else {
    names(ys) = NULL
  }

  return(ys)
}
