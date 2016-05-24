#' Generate or solve the cubic polynomial for ecd
#' 
#' Generate or solve the polynomial from ecd. 
#' This is usually transient for solve.
#' Or it can be used for studying singular points.
#'
#' @param object An object of ecd class
#' @param x A vector of x dimension
#' @param solve Logical, solve the polynomial, default = TRUE.
#'
#' @return list of the polynomial object, or result of solve.
#'
#' @keywords cubic
#'
#' @export
#'
#' @examples
#' d <- ecd()
#' ecd.cubic(d)
#' ecd.cubic(d, 0)
### <======================================================================>
"ecd.cubic" <- function(object, x=0, solve=TRUE)
{
    N <- length(x)
    
    a <- object@alpha
    r <- object@gamma
    s <- object@sigma
    b <- object@beta
    mu <- object@mu
    
    f <- function(x) {
        xi <- (x-mu)/s
        ec_poly <- polynom::polynomial(c(-a+xi^2, r+b*xi, 0, 1))
        if (solve) {
            return(solve(ec_poly))
        } else {
            return(ec_poly)
        }
    }
    lapply(x, f)
}
### <---------------------------------------------------------------------->
