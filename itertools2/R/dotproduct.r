#' Computes the dot product of two iterable objects
#'
#' Returns the dot product of two numeric iterables of equal length
#'
#' @importFrom iterators iter
#' @export
#' @param vec1 the first
#' @param vec2 the second iterable object
#' @return the dot product of the iterators
#'
#' @examples
#' it <- iterators::iter(1:3)
#' it2 <- iterators::iter(4:6)
#' dotproduct(it, it2) # 32
#'
#' it <- iterators::iter(1:4)
#' it2 <- iterators::iter(7:10)
#' dotproduct(1:4, 7:10) # 90
#'
dotproduct <- function(vec1, vec2) {
  vec1 <- iterators::iter(vec1)
  vec2 <- iterators::iter(vec2)
  it <- imap(prod, vec1, vec2)
  sum(sapply(it, identity))
}
