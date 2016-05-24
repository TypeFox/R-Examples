#' Sparsity
#' 
#' Show the sparsity (as a count or proportion) of a matrix.  For
#' example, .99 sparsity means 99\% of the values are zero.
#' Similarly, a sparsity of 0 means the matrix is fully dense.
#' 
#' @details
#' The implementation is very efficient for dense matrices.  For
#' sparse triplet matrices, the count is trivial.
#' 
#' @param x
#' The matrix, stored as an ordinary R matrix or as a "simple
#' triplet matrix" (from the slam package).
#' @param proportion
#' Logical; should a proportion or a count be returned?
#' 
#' @return
#' The sparsity of the input matrix, as a proportion or a count.
#' 
#' @examples
#' ## Completely sparse matrix
#' x <- matrix(0, 10, 10)
#' coop::sparsity(x)
#' 
#' ## 15\% density / 85\% sparsity
#' x[sample(length(x), size=15)] <- 1
#' coop::sparsity(x)
#' 
#' @author Drew Schmidt
#' @export
sparsity <- function(x, proportion=TRUE) UseMethod("sparsity")



#' @export
sparsity.matrix <- function(x, proportion=TRUE)
{
  if (is.integer(x))
    count <- .Call(R_sparsity_int, x)
  else if (is.double(x))
    count <- .Call(R_sparsity_dbl, x, tol=1e-10)
  else
    stop("matrix 'x' must be numeric.")
  
  if (proportion)
    count / nrow(x) / ncol(x)
  else
    count
}



#' @export
sparsity.simple_triplet_matrix <- function(x, proportion=TRUE)
{
  if (proportion)
    1 - length(x$v) / nrow(x) / ncol(x)
  else
    nrow(x)*ncol(x) - length(x$v)
}



# Sparse matrix generator; used only for tests
# @param m,n Dimensions (rows, cols)
# @param prop Proportion of non-zeros.
dense_stored_sparse_mat <- function(m, n, prop)
{
  size <- prop*m*n
  x <- matrix(0, m, n)
  x[sample(m*n, size=size)] <- 10#rnorm(size)
  x
}
