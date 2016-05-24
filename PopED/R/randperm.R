#' Function written to match MATLAB's randperm function
#' 
#' A wrapper for the \code{\link{sample}} function.
#' 
#' @param num Either a vector of one or more elements from which to choose, or a positive integer.
#' 
#' @return See \code{\link{sample}} 
#' 
#' @family MATLAB
#' @example tests/testthat/examples_fcn_doc/examples_randperm.R
#' @export
#' @keywords internal
## Function written to match MATLAB function
## Author: Andrew Hooker

randperm <- function(num){
    return(sample(num))
}
