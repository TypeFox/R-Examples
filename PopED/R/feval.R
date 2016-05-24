#' MATLAB feval function
#' 
#' This is just a wrapper for the \code{\link{do.call}} function to behave like the feval function in MATLAB.
#' 
#' @param file.name A function or a string that is the name of a function. 
#' @param ... Arguments for the function.  Multiple arguments separated by a comma.
#' 
#' @return Output from the defined function.
#' 
#' @family MATLAB
#' @example tests/testthat/examples_fcn_doc/examples_feval.R
#' @export
#' @keywords internal
#' 
## Function written to match MATLAB function
## Author: Andrew Hooker

feval <- function(file.name,...){
    #func.name <- gsub("\\.R$","",file.name)
    do.call(file.name,list(...))
}
