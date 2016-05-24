#' Iterator that returns an object indefinitely
#'
#' Constructs an iterator that returns an object over and over again.
#'
#' Runs indefinitely unless the \code{times} argument is specified. Used as
#' argument to \code{\link[itertools2]{imap}} for invariant function
#' parameters. Also used with \code{\link[itertools2]{izip}} to create constant
#' fields in a tuple record.
#'
#' @export
#' @param object object to return indefinitely.
#' @param times the number of times \code{object} is returned. If \code{NULL}
#' (default), \code{object} is returned indefinitely.
#' @return iterator that returns \code{object}
#' 
#' @examples
#' it <- irepeat(42)
#' iterators::nextElem(it)
#' iterators::nextElem(it)
#' iterators::nextElem(it)
#' # Further calls to iterators::nextElem(it) will repeat 42
#' 
#' it2 <- irepeat(42, times=4)
#' iterators::nextElem(it2)
#' iterators::nextElem(it2)
#' iterators::nextElem(it2)
#' iterators::nextElem(it2)
#'
#' # The object can be a data.frame, matrix, etc
#' it3 <- irepeat(iris, times=4)
#' iterators::nextElem(it3)
#' iterators::nextElem(it3)
#' iterators::nextElem(it3)
#' iterators::nextElem(it3)
irepeat <- function(object, times=NULL) {
  if (!is.null(times)) {
    times <- as.numeric(times)
    if (length(times) != 1) {
      stop("'times' must be a numeric value of length 1")
    }
  }

  i <- 0
  nextElem <- function() {
    i <<- i + 1

    if (!is.null(times) && i > times) {
      stop("StopIteration", call.=FALSE)
    }
    object
  }

  it <- list(nextElem=nextElem)
  class(it) <- c("abstractiter", "iter")
  it
}

