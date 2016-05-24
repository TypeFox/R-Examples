#' Consumes an iterator and computes its length
#'
#' Counts the number of elements in an iterator. NOTE: The iterator is consumed
#' in the process.
#'
#' @importFrom iterators iter nextElem
#' @export
#' @param object an iterable object
#' @return the number of elements in the iterator
#'
#' @examples
#' ilength(1:5) == length(1:5)
#' 
#' it <- iterators::iter(1:5)
#' ilength(it) == length(1:5)
#'
#' it2 <- ichain(1:3, 4:5, 6)
#' ilength(it2)
#'
#' it3 <- ichain(1:3, levels(iris$Species))
#' ilength(it3)
#'
ilength <- function(object) {
  it <- iterators::iter(object)

  i <- 0
  repeat{
    next_elem <- try(iterators::nextElem(it), silent=TRUE)
    if (stop_iteration(next_elem)) {
      break
    }
    i <- i + 1
  }
  i
}
