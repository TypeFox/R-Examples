#' Calculate range statistic.
#'
#' A simple wrapper for \code{diff(range(x))}, so \code{max(x) - min(x)}.
#'
#' @param x [\code{numeric}]\cr
#'   The vector.
#' @param na.rm [\code{logical(1)}]\cr
#'   If \code{FALSE}, NA is returned if an NA is encountered in \code{x}.
#'   If \code{TRUE}, NAs are disregarded.
#'   Default is \code{FALSE}
#' @return [\code{numeric(1)}].
#' @export
rangeVal = function(x, na.rm = FALSE) {
  assertNumeric(x, min.len = 1L, any.missing = TRUE)
  assertFlag(na.rm)
  if (allMissing(x))
    return(NA_real_)
  diff(range(x, na.rm = na.rm))
}
