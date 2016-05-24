#' Iterator that traverses each given iterable in a roundrobin order
#'
#' Constructs an iterator that traverses each given iterable in a roundrobin
#' order. That is, the iterables are traversed in an alternating fashion such
#' that the each element is drawn from the next iterable. If an iterable has no
#' more available elements, it is skipped, and the next element is taken from
#' the next iterable having available elements.
#' 
#' @importFrom iterators nextElem iter
#' @export
#' @param ... multiple arguments to iterate through in roundrobin sequence
#' @return iterator that alternates through each argument in roundrobin sequence
#'
#' @examples
#' it <- iterators::iter(c("A", "B", "C"))
#' it2 <- iterators::iter("D")
#' it3 <- iterators::iter(c("E", "F"))
#' as.list(iroundrobin(it, it2, it3)) # A D E B F C
#'
#' it_rr <- iroundrobin(1:3, 4:5, 7:10)
#' as.list(it_rr) # 1 4 7 2 5 8 3 9 10
#'
iroundrobin <- function(...) {
  iter_list <- lapply(list(...), iterators::iter)
  num_iters <- length(iter_list)

  has_elems <- rep(TRUE, num_iters)
  it_cycle <- icycle(seq_len(num_iters))

  nextElement <- function() {
    repeat {
      if (!any(has_elems)) {
        stop("StopIteration", call.=FALSE)
      }
      which_iter <- iterators::nextElem(it_cycle)
      next_elem <- try(iterators::nextElem(iter_list[[which_iter]]), silent=TRUE)

      if (stop_iteration(next_elem)) {
        has_elems[which_iter] <<- FALSE
      } else {
        return(next_elem)
      }
    }
  }

  it <- list(nextElem=nextElement)
  class(it) <- c("abstractiter", "iter")
  it
}
