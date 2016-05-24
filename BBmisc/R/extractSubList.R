#' Extracts a named element from a list of lists.
#'
#' @param xs [\code{list}]\cr
#'   A list of named lists.
#' @param element [\code{character}]\cr
#'   Name of element(s) to extract from the list elements of \code{xs}.
#'   What happens is this: \code{x$el1$el2....}.
#' @param element.value [any]\cr
#'   If given, \code{\link{vapply}} is used and this argument is passed to \code{FUN.VALUE}.
#'   Note that even for repeated indexing (if length(element) > 1) you only
#'   pass one value here which refers to the data type of the final result.
#' @param simplify [\code{logical(1)} | character(1)]\cr
#'   If \code{FALSE} \code{\link{lapply}} is used, otherwise \code{\link{sapply}}.
#'   If \dQuote{cols}, we expect the elements to be vectors of the same length and they are
#'   arranged as the columns of the resulting matrix.
#'   If \dQuote{rows}, likewise, but rows of the resulting matrix.
#'   Default is \code{TRUE}.
#' @param use.names [\code{logical(1)}]\cr
#'   If \code{TRUE} and \code{xs} is named, the result is named as \code{xs},
#'   otherwise the result is unnamed.
#'   Default is \code{TRUE}.
#' @return [\code{list} | simplified \code{vector} | \code{matrix}]. See above.
#' @export
#' @examples
#' xs = list(list(a = 1, b = 2), list(a = 5, b = 7))
#' extractSubList(xs, "a")
#' extractSubList(xs, "a", simplify = FALSE)
extractSubList = function(xs, element, element.value, simplify = TRUE, use.names = TRUE) {
  assertList(xs)
  assert(checkFlag(simplify), checkChoice(simplify, c("cols", "rows")))
  assertFlag(use.names)

  # we save some time here if we only do the for loop in the complicated case
  # the whole function is still not very nice due to the loop
  # extractSubList should be C code anyway i guess....
  doindex = if (length(element) == 1L) {
    function(x) x[[element]]
  } else {
    function(x) {
      for (el in element)
        x = x[[el]]
      return(x)
    }
  }

  if (!missing(element.value)) {
    ys = vapply(xs, doindex, FUN.VALUE = element.value)
  } else if (isTRUE(simplify)) {
    ys = sapply(xs, doindex, USE.NAMES = use.names)
   } else {
    ys = lapply(xs, doindex)
    if (simplify == "rows")
      ys = asMatrixRows(ys)
    else if (simplify == "cols")
      ys = asMatrixCols(ys)
  }
  ns = names(xs)
  if (use.names && !is.null(ns)) {
    if (isTRUE(simplify))
      names(ys) = ns
    else if (simplify == "rows")
      rownames(ys) = ns
    else if (simplify == "cols")
      colnames(ys) = ns
  } else {
    if (simplify %in% c("rows", "rows"))
      dimnames(ys) = NULL
    else
      names(ys) = NULL
  }
  return(ys)
}
