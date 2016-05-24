#' Iterator that replicates elements of an iterable object
#'
#' Constructs an iterator that replicates the values of an \code{object}.
#'
#' This function is intended an iterable version of the standard
#' \code{\link[base]{rep}} function. However, as exception, the recycling
#' behavior of \code{\link[base]{rep}} is intentionally not implemented.
#'
#' @export
#' @param object object to return indefinitely.
#' @param times the number of times to repeat each element in \code{object}
#' @param length.out non-negative integer. The desired length of the iterator
#' @param each non-negative integer. Each element is repeated \code{each} times
#' @return iterator that returns \code{object}
#' 
#' @examples
#' it <- irep(1:3, 2)
#' unlist(as.list(it)) == rep(1:3, 2)
#' 
#' it2 <- irep(1:3, each=2)
#' unlist(as.list(it2)) == rep(1:3, each=2)
#'
#' it3 <- irep(1:3, each=2, length.out=4)
#' as.list(it3)
irep <- function(object, times=1, length.out=NULL, each=NULL) {
  if (!is.null(length.out)) {
    length.out <- as.numeric(length.out)
    if (length(length.out) != 1) {
      stop("'length.out' must be a numeric value of length 1")
    }
  }
  if (!is.null(each)) {
    each <- as.numeric(each)
    if (length(each) != 1) {
      stop("'each' must be a numeric value of length 1")
    }
  }

  it <- icycle(object, times=times)
  if (!is.null(each)) {
    it <- irep_each(it, each=each)
  }

  islice(it, end=length.out)
}

#' @export
#' @rdname irep
irep_len <- function(object, length.out=NULL) {
  irep(object, times=1, length.out=length.out)
}

irep_each <- function(object, each=1) {
  each <- as.integer(each)
  iter_obj <- iterators::iter(object)
  
  iter_repeat <- irepeat(iterators::nextElem(iter_obj), times=each)

  nextElem <- function() {
    next_elem <- try(iterators::nextElem(iter_repeat), silent=TRUE)
    if (stop_iteration(next_elem)) {
      iter_repeat <<- irepeat(iterators::nextElem(iter_obj), times=each)
      next_elem <- iterators::nextElem(iter_repeat)
    }
    next_elem
  }

  it <- list(nextElem=nextElem)
  class(it) <- c("abstractiter", "iter")
  it
}
