#' Iterator that applies a given function to several iterables concurrently.
#'
#' Constructs an iterator that computes the given function \code{f} using the
#' arguments from each of the iterables given in \code{...}.
#'
#' The iterator returned is exhausted when the shortest iterable in \code{...}
#' is exhausted. Note that \code{imap} does not recycle arguments as
#' \code{\link[base]{Map}} does.
#'
#' The primary difference between \code{istarmap} and
#' \code{\link[itertools2]{imap}} is that the former expects an iterable object
#' whose elements are already grouped together, while the latter case groups the
#' arguments together before applying the given function. The choice is a matter
#' of style and convenience.
#'
#' @importFrom iterators nextElem
#' @export
#' @param f a function
#' @param ... multiple arguments to iterate through in sequence
#' @return iterator that returns the values of \code{object} along with the
#' index of the object. 
#' 
#' @examples
#' pow <- function(x, y) {
#'   x^y
#' }
#' it <- imap(pow, c(2, 3, 10), c(5, 2, 3))
#' as.list(it)
#'
#' # Similar to the above, but because the second vector is exhausted after two
#' # calls to `nextElem`, the iterator is exhausted.
#' it2 <- imap(pow, c(2, 3, 10), c(5, 2))
#' as.list(it2)
#'
#' # Another similar example but with lists instead of vectors
#' it3 <- imap(pow, list(2, 3, 10), list(5, 2, 3))
#' iterators::nextElem(it3) # 32
#' iterators::nextElem(it3) # 9
#' iterators::nextElem(it3) # 1000
imap <- function(f, ...) {
  f <- match.fun(f)
  iter_obj <- izip(...)

  nextElem <- function() {
    do.call(f, iterators::nextElem(iter_obj))
  }

  it <- list(nextElem=nextElem)
  class(it) <- c("abstractiter", "iter")
  it
}
