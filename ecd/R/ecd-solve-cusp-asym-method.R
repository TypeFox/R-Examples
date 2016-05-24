#' Trigonometric solution for asymmetric cusp distribution
#' 
#' The simplified trigonometric solution for \eqn{x^2=-y^3-beta*x*y}
#'
#'
#' @param x Array of x dimension
#' @param beta the skew parameter
#'
#' @return Array of y
#'
#' @keywords analytic
#'
#' @export
#'
#' @examples
#' x <- seq(-100,100,by=0.1)
#' y <- ecd.solve_cusp_asym(x, beta=0.5)

### <======================================================================>
"ecd.solve_cusp_asym" <- function(x, beta)
{
    if (length(beta)!=1) {
        stop("Asym cusp requires beta to be length-one numeric!")
    }
    
    if (beta==0) {
        return(-abs(x)^(2/3)) # revert to std cusp
    }

    # handle beta < 0
    if (beta < 0) {
        y <- ecd.solve_cusp_asym(-x, -beta)
        return(y)
    }
    
    # handle beta > 0
    x0 <- -(4*beta^3)/27
    V <- abs(x/x0)^(1/2)
    W <- 2*abs(beta*x/3)^(1/2)

    ifelse( x>=0, {
            -W*sinh(1/3*asinh(V))
        }, 
        ifelse( x<x0, {
            A <- suppressWarnings(acosh(V))
            -W*cosh(1/3*A)
        }, {
            A <- suppressWarnings(acos(V))
            -W*cos(1/3*A)
        })
    )

}
### <---------------------------------------------------------------------->
