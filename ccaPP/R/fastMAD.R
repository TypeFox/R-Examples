# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' Fast implementation of the median absolute deviation
#' 
#' Compute the median absolute deviation with a fast C++ implementation.  By 
#' default, a multiplication factor is applied for consistency at the normal 
#' model.
#' 
#' @param x  a numeric vector.
#' @param constant  a numeric multiplication factor.  The default value yields 
#' consistency at the normal model.
#' 
#' @return A list with the following components:
#' @returnItem center  a numeric value giving the sample median.
#' @returnItem MAD  a numeric value giving the median absolute deviation.
#' 
#' @note Functionality for removing observations with missing values is 
#' currently not implemented.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{fastMedian}}, \code{\link[stats]{mad}}
#' 
#' @examples 
#' set.seed(1234)  # for reproducibility
#' x <- rnorm(100)
#' fastMAD(x)
#' 
#' @keywords multivariate robust
#' 
#' @importFrom Rcpp evalCpp
#' @useDynLib ccaPP
#' @export

fastMAD <- function(x, constant = 1.4826) {
    # initializations
    x <- as.numeric(x)
    if(length(x) == 0) return(NA)  # zero length vector
    constant <- as.numeric(constant)
    # call C++ function
    .Call("R_fastMAD", R_x=x, R_constant=constant, PACKAGE="ccaPP")
}
