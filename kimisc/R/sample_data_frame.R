#' @title Random Samples and Permutations for Data Frames
#' @description This function takes a sample of the specified size
#'   from the rows of \code{x} using either with or without replacement.
#' @param x A data frame.
#' @param size A non-negative integer giving the number of items to choose.
#' @param replace Should sampling be with replacement?
#' @param prob A vector of probability weights for obtaining the rows of the
#'   data frame being sampled.
#' @return A data frame of the same shape as \code{x}.
#' @details This function internally calls \code{\link{sample.int}}.
#' @examples
#' set.seed(42)
#' sample.rows(data.frame(a=c(1,2,3), b=c(4,5,6), row.names=c('a', 'b', 'c')), 10, replace=TRUE)
#' @export
#' @rdname sample-methos
sample.rows <- function(x, size, replace=FALSE, prob=NULL) {
  stopifnot(is.data.frame(x))
  x.indexes <- sample.int(nrow(x), size, replace, prob)
  x[x.indexes, ]
}
