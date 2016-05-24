#' Creates a matrix of ones
#' 
#' Function creates a matrix of ones of size (dim1 x dim2). Written to match MATLAB's \code{ones} function.
#' 
#' 
#' @param dim1 The dimension of the matrix (if square) or the number of rows.   
#' @param dim2 The number of columns 
#' @return A matrix of ones
#' @family MATLAB
#' @example tests/testthat/examples_fcn_doc/examples_ones.R
#' @export
## Function written to match MATLAB function
## Author: Andrew Hooker

ones <- function(dim1,dim2=NULL){
    if(is.null(dim2)) dim2 <- dim1
    mat <- matrix(1,dim1,dim2)
    return(mat)
}
