#' Convert a numerical vector into a range string.
#'
#' @param x [\code{integer}]\cr
#'   Vector to convert into a range string.
#' @param range.sep [\code{character(1)}]\cr
#'   Separator between the first and last element of a range of consecutive
#'   elements in \code{x}.
#'   Default is \dQuote{ - }.
#' @param block.sep [\code{character(1)}]\cr
#'   Separator between non consecutive elements of \code{x} or ranges.
#'   Default is \dQuote{, }.
#' @return [\code{character(1)}]
#' @examples
#' x = sample(1:10, 7)
#' toRangeStr(x)
#' @export
toRangeStr = function(x, range.sep = " - ", block.sep = ", ") {
  if (testIntegerish(x))
    x = as.integer(x)
  else
    assertNumeric(x, any.missing = FALSE)
  assertString(range.sep)
  assertString(block.sep)

  findRange = function(x) seq_len(max(which(x == x[1L] + 0:(length(x)-1L))))
  x = sort(unique(x))
  x = unname(split(x, c(0L, cumsum(diff(x) > 1L))))
  combine = function(x)
    if (length(x) == 1L)
      as.character(x)
    else
      sprintf("%i%s%i", x[1L], range.sep, x[length(x)])
  collapse(vapply(x, combine, character(1L), USE.NAMES = FALSE), block.sep)
}
