#' Iterator that returns the elements of an object along with their indices
#'
#' Constructs an iterator that returns the elements of an object along with each
#' element's indices. Enumeration is useful when looping through an
#' \code{object} and a counter is required.
#'
#' This function is intended to follow the convention used in Python's
#' \code{enumerate} function where the primary difference is that a list is
#' returned instead of Python's \code{tuple} construct.
#'
#' Each call to \code{\link[iterators]{nextElem}} returns a list with two
#' elements:
#' \describe{
#'   \item{index:}{a counter}
#'   \item{value:}{the current value of \code{object}}
#' }
#'
#' \code{ienum} is an alias to \code{ienumerate} to save a few keystrokes.
#'
#' @export
#' @param object object to return indefinitely.
#' @return iterator that returns the values of \code{object} along with the
#' index of the object. 
#' 
#' @examples
#' set.seed(42)
#' it <- ienumerate(rnorm(5))
#' as.list(it)
#'
#' # Iterates through the columns of the iris data.frame
#' it2 <- ienum(iris)
#' iterators::nextElem(it2)
#' iterators::nextElem(it2)
#' iterators::nextElem(it2)
#' iterators::nextElem(it2)
#' iterators::nextElem(it2)
ienumerate <- function(object) {
  izip(index=icount(start=1), value=object)
}

#' @rdname ienumerate
#' @export
ienum <- ienumerate
