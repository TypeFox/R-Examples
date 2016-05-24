## Function written to match MATLAB function
## Author: Andrew Hooker

#' Compute the inverse of a matrix
#'
#' Function computes the inverse of a matrix.
#'
#' @param mat A matrix
#'
#' @return The inverse matrix
#' @export
#' @keywords internal

inv<- function(mat){
    #return(solve(mat))
    return(chol2inv(chol(mat)))
}
