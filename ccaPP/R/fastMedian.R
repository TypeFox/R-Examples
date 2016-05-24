# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' Fast implementation of the median
#' 
#' Compute the sample median with a fast C++ implementation.
#' 
#' @param x  a numeric vector.
#' 
#' @return The sample median.
#' 
#' @note Functionality for removing observations with missing values is 
#' currently not implemented.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{fastMAD}}, \code{\link[stats]{median}}
#' 
#' @examples 
#' set.seed(1234)  # for reproducibility
#' x <- rnorm(100)
#' fastMedian(x)
#' 
#' @keywords multivariate robust
#' 
#' @importFrom Rcpp evalCpp
#' @useDynLib ccaPP
#' @export

fastMedian <- function(x) {
    # initializations
    x <- as.numeric(x)
    if(length(x) == 0) return(NA)  # zero length vector
    # call C++ function
    .Call("R_fastMedian", R_x=x, PACKAGE="ccaPP")
}
