#' Iterator that returns elements of an object in pairs
#' 
#' Constructs an iterator of an iterable \code{object} that returns its elements
#' in pairs.
#' 
#' @importFrom iterators nextElem iter
#' @export
#' @param object an iterable object
#' @return an iterator that returns pairwise elements
#' 
#' @examples
#' it <- ipairwise(iterators::iter(letters[1:4]))
#' iterators::nextElem(it) # list("a", "b")
#' iterators::nextElem(it) # list("b", "c")
#' iterators::nextElem(it) # list("c", "d")
#'
#' it2 <- ipairwise(1:5)
#' iterators::nextElem(it2) # list(1, 2)
#' iterators::nextElem(it2) # list(2, 3)
#' iterators::nextElem(it2) # list(3, 4)
#' iterators::nextElem(it2) # list(4, 5)
#'
ipairwise <- function(object) {
  it_tee <- itee(object, n=2)
  dev_null <- iterators::nextElem(it_tee[[2]])
  
  nextElement <- function() {
    lapply(it_tee, iterators::nextElem)
  }

  it <- list(nextElem=nextElement)
  class(it) <- c("abstractiter", "iter")
  it
}

