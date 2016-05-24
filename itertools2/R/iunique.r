#' Iterator that extracts the unique elements from an iterable object
#'
#' Constructs an iterator that extracts each unique element in turn from an
#' iterable \code{object}. Order of the elements is maintained. This function is
#' an iterator analogue to \code{\link[base]{sort}}.
#'
#' NOTE: In order to determine whether an element is unique, a list of previous
#' unique elements is stored. In doing so, the list can potentially become large
#' if there are a large number of unique elements.
#'
#' @importFrom iterators nextElem
#' @export
#' @param object an iterable object
#' @return an iterator that returns the unique elements from \code{object}
#'
#' @examples
#' it <- ichain(rep(1, 4), rep(2, 5), 4:7, 2)
#' as.list(iunique(it)) # 1 2 4 5 6 7
#'
#' it2 <- iterators::iter(c('a', 'a', "A", "V"))
#' as.list(iunique(it2)) # a A V
#'
#' x <- as.character(gl(5, 10))
#' it_unique <- iunique(x)
#' as.list(it_unique) # 1 2 3 4 5
#'
iunique <- function(object) {
  iter_object <- iterators::iter(object)
  unique_elems <- list()
  i <- 1

  nextElem <- function() {
    repeat {
      next_elem <- try(iterators::nextElem(iter_object), silent=TRUE)
      if (stop_iteration(next_elem)) {
        stop("StopIteration", call.=FALSE)
      }

      if (!(next_elem %in% unique_elems)) {
        unique_elems[[i]] <<- next_elem
        i <<- i + 1
        return(next_elem)
      }
    }
  }

  it <- list(nextElem=nextElem)
  class(it) <- c("abstractiter", "iter")
  it
}

#' Iterator that extracts the just-seen unique elements from an iterable object
#'
#' Constructs an iterator that extracts each unique element in turn from an
#' iterable \code{object}. Order of the elements is maintained. Only the element
#' just seen is remembered for determining uniqueness.
#'
#' @importFrom iterators nextElem iter
#' @export
#' @param object an iterable object
#' @return an iterator that returns the just-seen unique elements from
#' \code{object}
#'
#' @examples
#' it <- ichain(rep(1,4), rep(2, 5), 4:7, 2)
#' it_iunique <- iunique_justseen(it)
#' as.list(it_iunique) # 1 2 4 5 6 7 2
#'
#' it2 <- iterators::iter(c('a', 'a', "A", 'a', 'a', "V"))
#' it2_iunique <- iunique_justseen(it2)
#' as.list(it2_iunique) # a A a V
#'
iunique_justseen <- function(object) {
  iter_object <- iterators::iter(object)
  prev_elem <- NULL
  first_seen <- FALSE

  nextElem <- function() {
    repeat {
      next_elem <- iterators::nextElem(iter_object)
      if (!first_seen || next_elem != prev_elem) {
        first_seen <<- TRUE
        prev_elem <<- next_elem
        return(next_elem)
      }
      prev_elem <<- next_elem
    }
  }

  it <- list(nextElem=nextElem)
  class(it) <- c("abstractiter", "iter")
  it
}
