#' Return the first n elements of an iterable object as a list
#'
#' Returns the first \code{n} elements of an iterable \code{object} as a list.
#' If \code{n} is larger than the number of elements in \code{object}, the
#' entire iterator is consumed.
#'
#' @export
#' @param object an iterable object
#' @param n the number of elements to return in the list
#' @return a list of the first \code{n} items of the iterable \code{object}
#'
#' @examples
#' take(iterators::iter(1:10), 3) # 1 2 3
#'
#' take(iterators::iter(1:5), 10) # 1 2 3 4 5
#'
take <- function(object, n=1) {
  if (n < 0 | !is.numeric(n) | length(n) != 1) {
    stop("n must be a positive integer of length 1")
  }
  n <- as.integer(n)
  as.list(islice(object, end=n))
}
