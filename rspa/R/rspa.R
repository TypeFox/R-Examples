#' A package for minimal vector adjustment.
#' 
#' @section Overview:
#'
#' {Given a vector \eqn{\boldsymbol{x}^0}, and a set linear restrictions of the
#' form  \eqn{\boldsymbol{a}_i\cdot\boldsymbol{x}_i=b_i} and/or 
#' \eqn{\boldsymbol{a}_i\cdot\boldsymbol{x}_i\leq b_i} with \eqn{i=1,2,\ldots,m}. This package finds the 
#' nearest vector to \eqn{\boldsymbol{x}^0} (in the (weighted) euclidean sense)
#' that obeys all restrictions. 
#' 
#' The package can handle large (sparse) problemes and has been tested by us on
#' vectors with half a million variables and sixty-thousand constraints. All underlying
#' algorithms have been implemented as a \code{C} library.
#' }
#' @name rspa-package
#' @docType package
#' @useDynLib rspa
#' @import editrules graphics
#' @importFrom stats density
{}


