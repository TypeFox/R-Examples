#' Sort the rows of a data.frame according to one or more columns.
#'
#' @param x [\code{data.frame}]\cr
#'   Data.frame to sort.
#' @param col [\code{character}]\cr
#'   One or more column names to sort \code{x} by.
#'   In order of preference.
#' @param asc [\code{logical}]\cr
#'   Sort ascending (or descending)?
#'   One value per entry of \code{col}.
#'   If a scalar logical is passed, it is replicated.
#'   Default is \code{TRUE}.
#' @return [\code{data.frame}].
#' @export
sortByCol = function(x, col, asc = TRUE) {
  assertDataFrame(x)
  assertSubset(col, colnames(x))
  m = length(col)
  assertLogical(asc, min.len = 1L, any.missing = FALSE)
  if (length(asc) == 1L)
    asc = rep(asc, m)

  asc = ifelse(asc, 1, -1)
  args = as.list(x[, col, drop = FALSE])
  # convert col to orderable numeric and multiply with factor
  args = Map(function(a, b) xtfrm(a) * b, args, asc)
  # now order the numerics and permute df
  o = do.call(order, args)
  return(x[o, , drop = FALSE])
}
