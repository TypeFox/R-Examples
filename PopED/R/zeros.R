#' Creates a matrix of zeros.
#' 
#' Function creates a matrix of zeros of size (dim1 x dim2). Written to match MATLAB's \code{zeros} function.
#' 
#' 
#' @param dim1 The dimension of the matrix (if square) or the number of rows.   
#' @param dim2 The number of columns 
#' @return A matrix of zeros.
#' @family MATLAB
#' @example tests/testthat/examples_fcn_doc/examples_zeros.R
#' @export
## Function written to match MATLAB function
## Author: Andrew Hooker

zeros <- function(dim1,dim2=NULL){
    if(is.null(dim2)){
        if(length(dim1)==2){
            tmp <- dim1
            dim1 <- tmp[1]
            dim2 <- tmp[2]
        } else if(length(dim1)==1){
            dim2 <- dim1
        } else {
            stop("first argument can only have one or two values")
        }
    }

    mat <- matrix(0,dim1,dim2)
    return(mat)
}
