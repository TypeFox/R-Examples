#' Function written to match MATLAB's isfield function
#' 
#' Check if a list or dataframe has an element with a specific name.
#' 
#' @param obj A list or dataframe    
#' @param sub.obj.str A string giving the name of the sub-object you want to check for.  
#' @return Logical. True if the element exists.
#' 
#' @family MATLAB
#' @example tests/testthat/examples_fcn_doc/examples_isfield.R
#' @export
#' @keywords internal
## Function written to match MATLAB function
## Author: Andrew Hooker

isfield <- function(obj,sub.obj.str){
    !is.null(obj[[sub.obj.str]])
}
