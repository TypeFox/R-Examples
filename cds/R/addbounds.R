#' Augment with Boundaries Between Rating Scale Categories and Rank
#' 
#' Adds \code{q - 1} boundaries between the \code{q} ratings to the columns of
#' matrix \code{x}, and convert the rows to rankings, starting with 0 for the
#' lowest ranking. Ties are handled by averaging the total rank for all 
#' tied observations.
#' 
#' Any \code{x} which is not a matrix or data frame will cause an error.
#' 
#' @param x matrix (or data frame) of \code{n} rows and \code{m} columns, or an
#' object that can be coerced to a matrix via \code{\link{as.matrix}}.
#' @param q scalar; the number of rating scale categories. Defaults to the
#' maximum entry in \code{x}.
#' @param ties character; handling of ties in \code{rank}
#' @return A matrix of size \code{n} by \code{m + q - 1}
#' @author Pieter C. Schoonees
#' @examples
#' set.seed(1234)
#' mat <- matrix(sample(1:9, 12, replace = TRUE), nrow = 4, ncol = 3)
#' addbounds(mat, q = 9)
#' @export
addbounds <- function(x, q = max(x), ties = "average"){
  if(is.data.frame(x)) x <- as.matrix(x)
  if(!is.matrix(x))
    stop("x must be a matrix or data frame")
  n <- nrow(x)
  m <- ncol(x)
  boundmat <- matrix(1:(q-1) + 0.5, byrow = TRUE, ncol = q - 1, nrow = n)
  out <- t(apply(cbind(x, boundmat), 1, rank, ties = ties)) - 1
  if(!is.null(colnames(x))) colnames(out) <- c(colnames(x), paste0("b", 1:(q - 1)))
  else colnames(out) <- c(paste0("Item", 1:m), paste0("b", 1:(q - 1)))
  out
}
