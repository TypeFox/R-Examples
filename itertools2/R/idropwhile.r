#' Iterator that drops elements until the predicate function returns FALSE
#'
#' Constructs an iterator that drops elements from the iterable \code{object} as
#' long as the \code{predicate} function is true; afterwards, every element of
#' \code{iterable} object is returned.
#'
#' Because the iterator does not return any elements until the \code{predicate}
#' first becomes false, there may have a lengthy start-up time before elements
#' are returned.
#'
#' @importFrom iterators iter nextElem
#' @export
#' @param predicate a function that determines whether an element is \code{TRUE}
#' or \code{FALSE}. The function is assumed to take only one argument.
#' @param object an iterable object
#' @return iterator object
#' 
#' @examples
#' # Filters out numbers exceeding 3
#' not_too_large <- function(x) {
#'   x <= 3
#' }
#' it <- idropwhile(not_too_large, 1:8)
#' as.list(it)
#'
#' # Same approach but uses an anonymous function
#' it2 <- idropwhile(function(x) x <= 10, seq(2, 20, by=2))
#' as.list(it2)
idropwhile <- function(predicate, object) {
  iter_obj <- iterators::iter(object)

  nextElem <- function() {
    repeat {
      next_elem <- iterators::nextElem(iter_obj)
      if (!predicate(next_elem)) {
        return(next_elem)
      }
    }
  }

  it <- list(nextElem=nextElem)
  class(it) <- c("abstractiter", "iter")
  it
}
