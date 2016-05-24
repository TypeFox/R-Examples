# :vim set filetype=R
#' Sample sub-sequences from a sequence
#'
#' This is like the normal sample function but instead of a scalar,
#' vector sub-sequences are extracted from the input.
#'
#' @section Usage:
#' samplerange(x, size, window, ...)
#'
#' @section Details:
#' Sometimes a sequence is auto-correlated. Attempting to construct a
#' a sub-sequence by sampling from such a sequence will lose the 
#' auto-correlation embedded within the original sequence. The solution
#' is to draw random sub-sequences from the original sequence, which
#' is what this function does.
#'
#' This operation can be for both a sequence (i.e. a vector or array)
#' or a matrix/data.frame. If the latter, a sub-matrix is selected 
#' such that the columns of the matrix are preserved.
#' This behavior is consistent
#' with time series data formats where a single series is represented
#' by a column and each row represents a point in time. Hence, the
#' 2D version will select sub-sequences in time, collecting all
#' associated time series.
#'
#' Under the hood, this function relies on sample.int, so the behavior
#' of the output can be controlled by passing additional arguments to
#' sample.int, such as replace=TRUE.
#'
#' @name samplerange
#' @param x A one-dimensional or two-dimensional data structure
#' @param size The number of sub-sequences to create
#' @param window The length of the output vectors 
#' @param \dots Optional arguments for the sample.int function
#'
#' @return 
#' When a sequence is passed to samplerange a matrix is returned,
#' where each column represents a sampled subsequence. Hence the
#' dimensions of the matrix will be window by size.
#'
#' If a matrix is passed to samplerange then a list of sub-matrices
#' is returned. Each sub-matrix will be of dimension window by ncol(x).
#' The length of the resulting list will be size.
#'
#' In either case, each _column_ is independent.
#'
#' @author Brian Lee Yung Rowe
#'
#' @examples
#' # Extract seven sub-sequences, each with length 3
#' samplerange(1:20, 7, 3)
#'
#' # This time use replacement
#' samplerange(1:20, 7, 3, replace=TRUE)
#'
#' # Extract five sub-matrices with dimensions 2 by 4
#' samplerange(matrix(1:32, ncol=4), 5, 2)
samplerange(x, size, window, ...) %when% {
  is.null(dim(x))
  window < length(x)
} %as% {
  count <- length(x) - window + 1
  samples <- sample.int(count, size, ...)
  sapply(samples, function(s) x[s:(s+window-1)])
}

samplerange(x, size, window, ...) %when% {
  window < length(x)
} %as% {
  count <- nrow(x) - window + 1
  samples <- sample.int(count, size, ...)
  lapply(samples, function(s) x[s:(s+window-1),])
}

