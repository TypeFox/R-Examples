#' Trigonometric solution for a elliptic curve
#' 
#' Use Chebyshev trigonometry for a depressed cube to solve a elliptic curve \eqn{y(x)}.
#'
#' @method solve_trig ecd
#'
#' @param object an object of ecd class
#' @param x array of x dimension
#'
#' @return array of y
#'
#' @keywords solve
#'
#' @author Stephen H-T. Lihn
#'
#' @export solve_trig
#'
#' @examples
#' d <- ecd()
#' x <- seq(-100,100,by=0.1)
#' y <- solve_trig(d,x)

### <======================================================================>
"solve_trig.ecd" <- function(object, x)
{
    cusp <- object@cusp
        
    xi <- (x-object@mu)/object@sigma
    alpha <- object@alpha - xi^2
    gamma <- object@gamma + object@beta*xi
    
    discr <- -16*(4*gamma^3 + 27*alpha^2)
    discr <- ifelse(cusp > 0 & abs(discr) < 1000*.Machine$double.eps, 0, discr)
    
    sqrt_3g <- sqrt(3/abs(gamma))/gamma
    sqrt_g3 <- sqrt(abs(gamma)/3)
    
    ifelse (gamma == 0, {
        # j=0 line
        sign(alpha)*abs(alpha)^(1/3)
    }, ifelse (discr >= 0, {
        # lower central region, also implies gamma <= 0
        # -C(-d@alpha, d@gamma)
        V <- 3/2*alpha*sqrt_3g
        A <- suppressWarnings(acos(V))
        -2*sqrt_g3 * cos(1/3*A)
        
    }, ifelse (gamma < 0 & discr < 0, {
        # Discr < 0, the outer lower left region
        # gamma < 0
        V <- -3/2*abs(alpha)*sqrt_3g
        A <- suppressWarnings(acosh(V))
        2*sign(alpha)*sqrt_g3 * cosh(1/3*A)
        
    }, { # discr < 0, also implies gamma >= 0
        # entire upper plane
        # but when gamma == 0, sqrt_3g and sqrt_g3 crap out
        V <- -3/2*alpha*sqrt_3g
        A <- suppressWarnings(asinh(V))
        -2*sqrt_g3 * sinh(1/3*A)
    })))
    
}
### <---------------------------------------------------------------------->
#' @rdname solve_trig.ecd
setGeneric("solve_trig", function(object, x) standardGeneric("solve_trig"))
#' @rdname solve_trig.ecd
setMethod("solve_trig", signature("ecd"), solve_trig.ecd)
### <---------------------------------------------------------------------->
