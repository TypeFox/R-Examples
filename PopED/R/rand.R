#' Function written to match MATLAB's rand function
#' 
#' Generate random samples from a uniform distribution [0,1] and return in matrix form
#' 
#' @param dim1 The dimension of the matrix (if square), otherwise the number of rows.
#' @param dim2 The number of colums, if different from the number of rows.
#' 
#' @return Matrix of random generated samples.
#' 
#' @family MATLAB
#' @example tests/testthat/examples_fcn_doc/examples_rand.R
#' @export
#' @keywords internal

## Function written to match MATLAB function
## Author: Andrew Hooker

rand <- function(dim1,dim2=NULL){
    if(is.null(dim2)) dim2 <- dim1
    tmp <- runif(dim1*dim2)
    mat <- matrix(tmp,dim1,dim2)
    return(mat)
}
