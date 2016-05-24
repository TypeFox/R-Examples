#' Bootstrap Resampling
#'
#' Generate a dataset from a given dataset using sampling with or without replacement.
#'
#' This function is a simple shortcut for generating a bootstrap sample. It was
#' originally designed for use within the comparison framework, but can be used as a
#' stand-alone function.
#' @param dataset a p x m data matrix.
#' @param n the number of rows of observations to be included in the new dataset.
#' @param replace logical. If replace is TRUE sampling will be with replacement, thus
#'        returning a bootstrap replicate.
#' @return The function returns a new n x m matrix. The default is a bootstrap
#'         replicate of the orginal data.

randboot <- function(dataset, n, replace = TRUE) {

  m <- dim(dataset)[1]

  if (missing(n)) n <-  m
  y <- dataset[sample(m, n, replace), ]

  return(y)

}
