#' Iterator that returns elements in fixed-length chunks
#'
#' Constructs an iterator that returns elements of an iterable \code{object} in
#' fixed-length chunks. If the length of the iterator is not divisible by
#' \code{chunk_size}, the remainder of the last block is filled with the value
#' specified in \code{fill}.
#'
#' This function corresponds to Python's \code{grouper} function. We chose the
#' name \code{ichunk} because it more explicitly defines the function's purpose.
#'
#' @importFrom iterators nextElem iter
#' @export
#' @param object an iterable object
#' @param chunk_size the number of elements returned per chunk
#' @param fill the value with which to fill the last chunk if the length of the
#' iterator is not divisble by \code{chunk_size}
#' @return each call to \code{nextElem} results in a list of length
#' \code{chunk_size}
#'
#' @examples
#' it <- ichunk(iterators::iter(1:5), chunk_size=2)
#' # List: list(1, 2, 3)
#' iterators::nextElem(it)
#' # List: list(4, 5, NA)
#' iterators::nextElem(it)
#'
#' it2 <- ichunk(levels(iris$Species), chunk_size=4, "weeee")
#' # Returns: list("setosa", "versicolor", "virginica", "weeee")
#' iterators::nextElem(it2)
#'
ichunk <- function(object, chunk_size=1, fill=NA) {
  if (chunk_size <= 0 | !is.numeric(chunk_size) | length(chunk_size) != 1) {
    stop("'chunk_size' must be a positive integer of length 1")
  }
  chunk_size <- as.integer(chunk_size)

  it <- iterators::iter(object)
  it_replicate <- replicate(n=chunk_size, it, simplify=FALSE)

  nextElem <- function() {
    Map(try_nextElem, it_replicate, default=fill)
  }

  it_chunk <- list(nextElem=nextElem)
  class(it_chunk) <- c("abstractiter", "iter")
  it_chunk
}
