#' @title Create Sparse Matrix
#' 
#' @description Creates a sparse square matrix with a given sparsity and distribution.
#' 
#' @param N the dimension of the square matrix
#' @param sparsity the density of non zero elements
#' @param method the method used to generate the entries of the matrix. Possible values are 
#' \code{"normal"} (default) or \code{"bimodal"}.
#' @param stationary should the spectral radius of the matrix be smaller than 1? 
#' Possible values are \code{TRUE} or \code{FALSE}. Default is \code{FALSE}.
#' @return An NxN sparse matrix. 
#' @examples
#' M <- createSparseMatrix(N = 30, sparsity = 0.05, method = "normal", stationary = TRUE)
#'
#' @export
createSparseMatrix <- function(N, sparsity, method = "normal", stationary = FALSE) {
  
  if (method == "normal") {
    # normal distributed nonzero entries
    n <- floor(sparsity * N^2)
    nonZeroEntries <- stats::rnorm(n, mean = 0, sd = 1)
    entries <- sample(x = 1:N^2, size = n, replace = FALSE)
    
    Atmp <- numeric(length = N^2)
    Atmp[entries] <- nonZeroEntries
    
    A <- matrix(Atmp, nrow = N, ncol = N)
    
  } else if (method == "unimodal") {
    
    # normal distributed nonzero entries
    n <- floor(sparsity * N^2)
    nonZeroEntries <- stats::rnorm(n, mean = 0, sd = 1)
    entries <- sample(x = 1:N^2, size = n, replace = FALSE)
    
    Atmp <- numeric(length = N^2)
    Atmp[entries] <- nonZeroEntries
    
    A <- matrix(Atmp, nrow = N, ncol = N)
    
  } else if (method == "bimodal") {
    # bimodal distributed nonzero entries
    n <- floor(sparsity * N^2)
    
    nonZeroEntriesLeft <- stats::rnorm(n, mean = -1, sd = sqrt(0.5))
    nonZeroEntriesRight <- stats::rnorm(n, mean = 1, sd = sqrt(0.5))
    
    nonZeroEntries <- sample(x = c(nonZeroEntriesLeft, nonZeroEntriesRight), size = n, replace = FALSE)
    
    entries <- sample(x = 1:N^2, size = n, replace = FALSE)
    
    Atmp <- numeric(length = N^2)
    Atmp[entries] <- nonZeroEntries
    
    A <- matrix(Atmp, nrow = N, ncol = N)

  } else {
    # invalid method
    stop("Unknown method. Possible methods are normal or bimodal.")
    
  }

  if (stationary == TRUE) {
    # if spectral radius < 1 is needed, return the re-normalize the matrix  
    K <- 1
    return(1/(K * sqrt(sparsity * N)) * A)
    #return(1/(max(Mod(eigen(A)$values)) + 0.01) * A)
    
  } else {
    
    return(A)
    
  }
  
}
