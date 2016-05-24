#' Function written to match MATLAB's size function
#' 
#' @param obj An object you want to know the various dimensions of.  Typically a matrix.
#' @param dimension.index Which dimension you are interested in.
#' @return The dimensions of the object or specific dimension you are interested in. 
#' @family MATLAB
#' @example tests/testthat/examples_fcn_doc/examples_size.R
#' @export
## Function written to match MATLAB function
## Author: Andrew Hooker

size <- function(obj,dimension.index=NULL){
    dim.obj <- dim(obj)
    if(is.null(dim.obj)) dim.obj <- c(1,length(obj))
    ##if(is.null(dim.obj)) dim.obj <- c(length(obj),length(obj))
    if(is.null(dimension.index)) return(dim.obj)
    return(dim.obj[[dimension.index]])
}
