#' Returns a list of n independent iterators from a single iterable object
#'
#' Constructs a list of \code{n} iterators, each of which iterates through an
#' iterable \code{object}.
#'
#' If the \code{object} is an iterator (i.e., inherits from class \code{iter}),
#' \code{n} deep copies of \code{object} are returned. Otherwise, \code{object}
#' is passed to \code{\link[iterators]{iter}} \code{n} times.
#'
#' @export
#' @param object an iterable object
#' @param n the number of iterables to return
#' @return a list of \code{n} iterators
#' 
#' @examples
#' # Creates a list of three iterators.
#' # Each iterator iterates through 1:5 independently.
#' iter_list <- itee(1:5, n=3)
#'
#' # Consumes the first iterator
#' unlist(as.list(iter_list[[1]])) == 1:5
#'
#' # We can iterate through the remaining two iterators in any order.
#' iterators::nextElem(iter_list[[2]]) # 1
#' iterators::nextElem(iter_list[[2]]) # 2
#' 
#' iterators::nextElem(iter_list[[3]]) # 1
#' iterators::nextElem(iter_list[[3]]) # 2
#'
#' iterators::nextElem(iter_list[[2]]) # 3
#' iterators::nextElem(iter_list[[2]]) # 4
#' iterators::nextElem(iter_list[[2]]) # 5
#'
#' iterators::nextElem(iter_list[[3]]) # 3
#' iterators::nextElem(iter_list[[3]]) # 4
#' iterators::nextElem(iter_list[[3]]) # 5
itee <- function(object, n=2) {
  n <- as.integer(n)
  if (length(n) != 1) {
    stop("'n' must be an integer value of length 1")
  } else if (n < 1) {
    stop("'n' must be a positive integer")
  }

  # TODO: Confirm that efficient buffering is performed for a more complex
  # example. For assistance on this, see this excellent post about how Python's
  # itertools.tee works: http://discontinuously.com/2012/06/inside-python-tee/

  # If the 'object' is an iterator, n deep copies of 'object' are returned.
  # Otherwise, 'object' is passed to iterators::iter 'n' times.
  if (inherits(object, "iter")) {
    itee_list <- replicate(n=n,
                           expr=iter_deepcopy(object),
                           simplify=FALSE)
  } else {
    itee_list <- replicate(n=n,
                           expr=iterators::iter(object),
                           simplify=FALSE)
  }

  itee_list
}

