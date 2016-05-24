#' Consumes the first n elements of an iterator
#'
#' Advances the iterator n-steps ahead without returning anything.
#'
#' If \code{n} is 0, the iterator is consumed entirely. Similarly, if \code{n}
#' is larger than the length of the iterator, the iterator is consumed entirely.
#'
#' @importFrom iterators nextElem
#' @export
#' @param iterator an iterator object
#' @param n The number of elements to consume.
#' @return Nothing, i.e., \code{invisible(NULL)}
#'
#' @examples
#' it <- iterators::iter(1:10)
#' # Skips the first 5 elements
#' consume(it, n=5)
#' # Returns 6
#' iterators::nextElem(it)
#'
#' it2 <- iterators::iter(letters)
#' # Skips the first 4 elements
#' consume(it2, 4)
#' # Returns 'e'
#' iterators::nextElem(it2)
#' 
consume <- function(iterator, n=0) {
  if (!is_iterator(iterator)) {
    stop("'iterator' must be of class 'iter'")
  }
  
  if (n < 0 || !is.numeric(n) || length(n) != 1) {
    stop("n must be a non-negative integer of length 1")
  }

  n <- as.integer(n)
  i <- 0
  repeat {
    elem <- try(iterators::nextElem(iterator), silent=TRUE)
    i <- i + 1
    if (stop_iteration(elem) || (n > 0 && i >= n)) {
      break
    }
  }

  invisible(NULL)
}

#' Returns the nth item of an iterator
#'
#' Returns the \code{n}th item of an \code{iterator} after advancing the
#' iterator \code{n} steps ahead. If the \code{iterator} is entirely consumed,
#' the \code{default} value is returned instead. That is, if either \code{n >
#' length(iterator)} or \code{n} is 0, then the \code{iterator} is consumed.
#'
#' @importFrom iterators nextElem
#' @export
#' @param iterator an iterator object
#' @param n The location of the desired element to return
#' @param default The value to return if iterable is consumed, default is NA
#' @return The nth element of the iterable or the default value
#'
#' @examples
#' it <- iterators::iter(1:10)
#' # Returns 5
#' nth(it, 5)
#'
#' it2 <- iterators::iter(letters)
#' # Returns 'e'
#' nth(it2, 5)
#'
#' it3 <- iterators::iter(letters)
#' # Returns default value of NA
#' nth(it3, 42)
#'
#' it4 <- iterators::iter(letters)
#' # Returns default value of "foo"
#' nth(it4, 42, default="foo")
#'
nth <- function(iterator, n, default=NA) {
  if (!is_iterator(iterator)) {
    stop("'iterator' must be of class 'iter'")
  }

  if (n < 0 | !is.numeric(n) | length(n) != 1) {
    stop("n must be a positive integer of length 1")
  }

  n <- as.integer(n)

  it <- islice(iterator, start=n)
  elem <- try(iterators::nextElem(it), silent=TRUE)
  ifelse(stop_iteration(elem), default, elem)
}
