#' Iterator that iterates through several iterables concurrently.
#'
#' The resulting iterator aggregates elements from each of the iterables into a
#' list from each iteration. Used for lock-step iteration over several iterables
#' at a time.
#'
#' @importFrom iterators iter nextElem
#' @export
#' @param ... multiple arguments to iterate through in sequence
#' @return iterator that iterates through each argument in sequence
#' 
#' @examples
#' it <- izip(x=1:3, y=4:6, z=7:9)
#' iterators::nextElem(it) # list(x=1, y=4, z=7)
#' iterators::nextElem(it) # list(x=2, y=5, z=8)
#' iterators::nextElem(it) # list(x=3, y=6, z=9)
#'
#' # Sums the zip'd elements. 1 + 4 + 7, and so on.
#' it2 <- izip(1:3, 4:6, 7:9)
#' sum_zip <- sapply(it2, function(x) sum(unlist(x)))
#' sum_zip == c(12, 15, 18)
#' 
#' it3 <- izip(a=1:3, b=4:42, class=levels(iris$Species))
#' iterators::nextElem(it3) # list(a=1, b=4, class="setosa")
#' iterators::nextElem(it3) # list(a=2, b=5, class="versicolor")
#' iterators::nextElem(it3) # list(a=3, b=6, class="virginica")
izip <- function(...) {
  iter_list <- lapply(list(...), iter)
  if (length(iter_list) == 0) {
    stop("At least one argument must be supplied.")
  }

  # TODO: Introduce iterators:::hasNext() to circumvent catching StopIteration?
  #       See the iterators vignette for details.
  nextElem <- function() {
    next_elem <- lapply(iter_list, function(x) {
      try(iterators::nextElem(x), silent=TRUE)
    })

    any_stopped <- any(sapply(next_elem, stop_iteration))
    if (any_stopped) {
      stop("StopIteration", call.=FALSE)
    }
    next_elem
  }

  it <- list(nextElem=nextElem)
  class(it) <- c("abstractiter", "iter")
  it
}

