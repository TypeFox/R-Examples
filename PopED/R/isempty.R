#' Function written to match MATLAB's isempty function
#' 
#' @param ... arguments to pass to the function. Typically a matrix.   
#' @return Logical. True if the passed object has any dimension that is zero.
#' @family MATLAB
#' @example tests/testthat/examples_fcn_doc/examples_isempty.R
#' @export
#' @keywords internal
## Function written to match MATLAB function
## Author: Andrew Hooker

isempty <- function(...){
    ret <- any(size(...)==0)
    return(ret)
}
