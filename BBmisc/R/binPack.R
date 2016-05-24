#' Simple bin packing.
#'
#' Maps numeric items in \code{x} into groups with sum
#' less or equal than \code{capacity}.
#' A very simple greedy algorithm is used, which is not really optimized
#' for speed. This is a convenience function for smaller vectors, not
#' a competetive solver for the real binbacking problem.
#' If an element of \code{x} exceeds \code{capacity}, an error
#' is thrown.
#'
#' @param x [\code{numeric}]\cr
#'   Numeric vector of elements to group.
#' @param capacity [\code{numeric(1)}]\cr
#'   Maximum capacity of each bin, i.e., elements will be grouped
#'   so their sum does not exceed this limit.
#' @return [\code{integer}]. Integer with values \dQuote{1} to \dQuote{n.bins}
#'   indicating bin membership.
#' @export
#' @examples
#' x = 1:10
#' bp = binPack(x, 11)
#' xs = split(x, bp)
#' print(xs)
#' print(sapply(xs, sum))
binPack = function(x, capacity) {
  assertNumeric(x, min.len = 1L, lower = 0, any.missing = FALSE)
  assertNumber(capacity)

  too.big = which.first(x > capacity, use.names = FALSE)
  if (length(too.big))
    stopf("Capacity not sufficient. Item %i (x=%f) does not fit", too.big, x[too.big])
  if (any(is.infinite(x)))
    stop("Infinite elements found in 'x'")

  ord = order(x, decreasing = TRUE)
  grp = integer(length(x))
  sums = vector(typeof(x), 1L)
  bin.count = 1L

  for(j in ord) {
    new.sums = sums + x[j]
    pos = which.first(new.sums <= capacity, use.names = FALSE)
    if (length(pos)) {
      grp[j] = pos
      sums[pos] = new.sums[pos]
    } else {
      bin.count = bin.count + 1L
      grp[j] = bin.count
      sums[bin.count] = x[j]
    }
  }
  grp
}
