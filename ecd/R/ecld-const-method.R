#' Analytic solution of the normalization constant for lambda distribution
#' 
#' The normalization constant \eqn{C}. SGED is supported.
#'
#' @param object an object of ecld class
#'
#' @return numeric 
#'
#' @keywords ecld
#'
#' @author Stephen H. Lihn
#'
#' @export 
#'
#' @examples
#' ld <- ecld(3)
#' ecld.const(ld)

### <======================================================================>
"ecld.const" <- function(object)
{
    ecld.validate(object, sged.allowed=TRUE)
    one <- if(object@use.mpfr) ecd.mp1 else 1 # for gamma function

    lambda <- object@lambda
    s <- object@sigma
    b <- object@beta
    
    # SGED
    if (object@is.sged) {
        return(s*lambda*gamma(lambda/2))
    }

	# symmetric
    if (b==0) {
        return( s*lambda*gamma(lambda/2) )
    }

	# asymmetric
    if (lambda==2) {
        return( 2*s*sqrt(1+b^2/4) )
    }
    if (lambda > 2) {
        ld0 <- ecld(lambda=lambda, beta=b)
        d0 <- ecd(lambda=ecd.mp2f(lambda), beta=ecd.mp2f(b), sigma=one, bare.bone=TRUE)
        e_y2 <- function(x) {
            exp(ecld.solve(ld0,x)) + exp(ecld.solve(ld0,-x))
        }
        C <- ecd.integrate(d0, e_y2, 0, Inf)
        if (C$message != "OK") {
            stop("Failed to integrate C from unit distribution")
        }
        return(ecd.mpnum(object, s*C$value))
    }

    stop("Unknown analytic solution for const")
}
### <---------------------------------------------------------------------->

