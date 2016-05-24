#' Analytic solution for a symmetric elliptic curve
#' 
#' Analytic solution for a symmetric elliptic curve \eqn{y(x)}
#'
#' @method solve_sym ecd
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
#' @export solve_sym
#'
#' @examples
#' d <- ecd()
#' x <- seq(-100,100,by=0.01)
#' y <- solve_sym(d,x)

### <======================================================================>
"solve_sym.ecd" <- function(object, x)
{
    if(object@beta != 0) {
        stop("solve_sym.ecd cannot handle non-zero beta")
    }
    
    a <- object@alpha
    r <- object@gamma
    s <- object@sigma
    mu <- object@mu
    cusp <- object@cusp
    
    # ----------------
    xi <- (x-mu)/s

    if (r==0) { # this includes std cusp
        a1 <- a-xi^2
        return(sign(a1)*abs(a1)^(1/3))
    }    
    # ----------------
    f1 <- (xi^2 - a)/2 
    
    f2_core <- 27*xi^4 - 54*a*xi^2 + 27*a^2 + 4*r^3
    c2 <- 2*3^(3/2) 

    fsgn <- ifelse(f2_core >= 0, ifelse((f1*c2)^2 < f2_core, -1, 1), 1)
    
    f2 <- sqrt(f2_core+0i)/c2        
    f3 <- (fsgn*(f1-f2)+0i)^(1/3)
    
    # z0 <- (-1+0i)^(1/3)
    # z1 <- -0.5 + 0.5i*sqrt(3)
    # z2 <- -0.5 - 0.5i*sqrt(3) 
    # f4 <- z1*f3*z0 - r*z2/(3*f3*z0)
    # z1*z0 and z2/z0 are -1
    
    f5_complex <- (-f3 + r/(3*f3))*fsgn
    f5_pi <- Arg(f5_complex)/pi
    f5_is_real <- abs(f5_pi-round(f5_pi)) <= 0.0001
    f5 <- ifelse(f5_is_real, Re(f5_complex), NaN)
    
    f5
}
### <---------------------------------------------------------------------->
#' @rdname solve_sym.ecd
setGeneric("solve_sym", function(object, x) standardGeneric("solve_sym"))
#' @rdname solve_sym.ecd
setMethod("solve_sym", signature("ecd"), solve_sym.ecd)
### <---------------------------------------------------------------------->
