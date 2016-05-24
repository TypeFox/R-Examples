#' Analytic solution for \eqn{y(x)} in lambda distribution
#' 
#' Analytic solution for \eqn{y(x)} if available. 
#' \eqn{ecld.laplace_B} is a utility function for the slopes of 
#' a skew Laplace distribution at lambda=2: \eqn{B^+} and \eqn{B^-} 
#' with \eqn{B^0/2=B^{+}+B^{-}}. If \code{sigma} is provided, B notation
#' is expanded for IMGF where \eqn{B^{+}_{\sigma} B^{-}_{\sigma} = \exp(\mu_D)}.
#' SGED is supported.
#'
#' @param a an object of ecld class
#' @param b a vector of \eqn{x} values
#' @param sgn sign of \eqn{-1, 0, +1}
#' @param beta the skew parameter
#' @param sigma the scale parameter, optional
#' @param ... Not used. Only here to match the generic signature.
#'
#' @return A vector for \eqn{y(x)}
#'
#' @keywords solve
#'
#' @author Stephen H. Lihn
#'
#' @export ecld.solve
#' @export ecld.laplace_B
#'
#' @examples
#' ld <- ecld(sigma=0.01*ecd.mp1)
#' x <- seq(-0.1, 0.1, by=0.01)
#' ecld.solve(ld,x)

### <======================================================================>
"ecld.solve" <- function(a, b, ...)
{
	# convert to typical variable names
	object <- a
	x <- b

    ecld.validate(object, sged.allowed=TRUE)
    one <- if(object@use.mpfr) ecd.mp1 else 1 # for consistent type
    
    lambda <- object@lambda * one
    s <- object@sigma * one
    beta <- object@beta 
	xi <- (x-object@mu)/s
    
    # SGED
    if (object@is.sged) {
        s2 <- ecld.ifelse(object, xi<0, s*(1-beta), s*(1+beta))
        xi2 <- (x-object@mu)/s2
        return(-abs(xi2)^(2/lambda))
    }
	
	# symmetric
	if (beta==0) {
		return(-abs(xi)^(2/lambda))
	}
	
	# asymmetric
    if (lambda==2 & beta != 0) {
    	Bi <- ecld.laplace_B(beta, ifelse(xi<0, 1, -1))
        return(-Bi*abs(xi))
    }
    if (lambda==3 & beta != 0) {
    	x0 <- abs(4*beta^3/27)
    	V <- sqrt(abs(xi/x0))
    	W <- 2*sqrt(abs(xi*beta/3))
    	
    	I <- function(x) {
    	    if(beta > 0) sinh(x) else cosh(x)
    	}
    	J <- function(x) {
    	    if(beta > 0) cosh(x) else sinh(x)
    	}
    	Iarc <- function(x) {
    	    suppressWarnings( if(beta > 0) asinh(x) else acosh(x) )
    	}
    	Jarc <- function(x) {
    	    suppressWarnings( if(beta > 0) acosh(x) else asinh(x) )
    	}
    	
    	K <- function(x) {
    	    zero <- 0
    	    c <- ifelse(xi==0, zero, ifelse(xi>=0, I(x), J(x)))
    	    ecd.mpnum(object, Re(c))
    	}
    	Karc <- function(x) {
            xc <- x + 0i
    	    a1 <- ifelse(xi>=0, Iarc(x), Jarc(x))
    	    a2 <- ifelse(xi>=0, Iarc(xc), Jarc(xc))
            ifelse(is.na(a1), a2, a1)
    	}
        
        V1 <- ecd.mp2f(V) # can not use MPFR and complex together!
    	y <- -K(1/3*Karc(V1))*W
        if (any(is.na(y))) {
            xi_na <- xi[is.na(y)]
            stop(paste("NaN found in lambda=3 beta=", beta, "xi=", xi[1]))
        }
        return(y)
    }
	if (lambda > 2 & lambda != 3 & beta != 0) {
	    # enroute solve.ecd of unit distribution
	    d0 <- ecd(lambda=ecd.mp2f(lambda), beta=ecd.mp2f(beta), bare.bone=TRUE)
	    return(solve(d0,xi))
	}
	
    stop("Unknown analytic formula for CDF")

}
### <---------------------------------------------------------------------->
#' @rdname ecld.solve
ecld.laplace_B <- function(beta, sgn=0, sigma=0) {
	
    stopifnot(all(sgn==0 | abs(sgn)==1))
    stopifnot(all(!is.na(sigma)))
    stopifnot(all(!is.na(beta)))
    
	B0 <- sqrt(1+beta^2/4)
	return(B0 + sgn*beta/2 + sgn*sigma)
}