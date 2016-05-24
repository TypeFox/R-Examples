#' Estimate the normalization constant for an ecd object
#'
#' This is an internal helper function for ecd constructor
#' Its main function is to estimate \code{const} using analytical formula,
#' without any dependency on statistics and numerical integration.
#'
#' @param object An object of ecd class
#'
#' @return numeric, estimated const
#'
#' @keywords constructor
#'
#' @export
#'
#' @examples
#' ecd.estimate_const(ecd(100,100, sigma=0.1, bare.bone=TRUE))
#'
### <======================================================================>
ecd.estimate_const <- function(object) {
    R <- object@R
    y0 <- solve(object,0)
    var3 <- 3/2*R^(2/3) - 9/2*R^(1/3) + 63/8
    C3 <- 2*sqrt(pi/2) *exp(y0) *sqrt(var3) *object@sigma
    return(C3)
}
### <---------------------------------------------------------------------->

