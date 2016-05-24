#' Iterator of neverending numeric sequence with initial value and step size
#'
#' Constructs an iterator that generates a neverending sequence of evenly spaced
#' values starting with \code{icount}. The step size is given by \code{step}.
#'
#' NOTE: Use a negative \code{step} size to generate decreasing sequences.
#'
#' Often used as an argument to \code{\link[itertools2]{imap}} to
#' generate consecutive data points.
#'
#' @export
#' @param start sequence's initial value
#' @param step sequence's step size
#' @return sequence's iterator
#' 
#' @examples
#' it <- icount()
#' iterators::nextElem(it)
#' iterators::nextElem(it)
#' iterators::nextElem(it)
#' 
#' it2 <- icount(start=5.5, step=1.5)
#' iterators::nextElem(it2)
#' iterators::nextElem(it2)
#' iterators::nextElem(it2)
icount <- function(start=0, step=1) {
  start <- as.numeric(start)
  step <- as.numeric(step)

  if (length(start) != 1) {
    stop("'start' must be a numeric value of length 1")
  }
  if (length(step) != 1) {
    stop("'step' must be a numeric value of length 1")
  }

  current_val <- start - step
  nextElem <- function() {
    current_val <<- current_val + step
    current_val
  }

  it <- list(nextElem=nextElem)
  class(it) <- c("abstractiter", "iter")
  it
}

