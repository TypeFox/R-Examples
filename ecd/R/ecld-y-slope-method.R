#' Analytic solution for the slope of \eqn{y(x)} in lambda distribution
#' 
#' Analytic solution for the slope of \eqn{y(x)} if available.
#' \eqn{ecld.y_slope_trunc} calculates the MGF truncation point where
#' \eqn{dy/dx+t=1}. SGED is supported.
#'
#' @param object an object of ecld class
#' @param x a vector of \eqn{x} values
#' @param t numeric, for MGF truncation
#'
#' @return numeric
#'
#' @keywords y_slope
#'
#' @author Stephen H-T. Lihn
#'
#' @export ecld.y_slope
#' @export ecld.y_slope_trunc
#'
#' @importFrom stats uniroot
#'
#' @examples
#' ld <- ecld(sigma=0.01*ecd.mp1)
#' x <- seq(-0.1, 0.1, by=0.01)
#' ecld.y_slope(ld,x)
#' ecld.y_slope_trunc(ld)

### <======================================================================>
"ecld.y_slope" <- function(object, x)
{

    ecld.validate(object, sged.allowed=TRUE)
    one <- if(object@use.mpfr) ecd.mp1 else 1 # for consistent type
    
    lambda <- object@lambda * one
    s <- object@sigma * one
    beta <- object@beta 
	xi <- (x-object@mu)/s
	
    # SGED
    if (object@is.sged) {
        sgn <- sign(xi)
        sgn[xi==0] <- NaN
        s2 <- ecld.ifelse(object, xi<0, s*(1-beta), s*(1+beta))
        xi2 <- (x-object@mu)/s2
		return(-sgn * abs(xi2)^(2/lambda-1) *2/(lambda*s2))
	}
    
	# symmetric
	if (beta==0) {
        sgn <- sign(xi)
        sgn[xi==0] <- NaN
		return(-sgn * abs(xi)^(2/lambda-1) *2/(lambda*s))
	}

    # general
    y <- ecld.solve(object, x)
    p <- 2*xi + beta*y
    q <- -lambda*s*(-y)^(lambda-1) - beta*s*xi
    return(ecd.mpnum(object,p/q))
	
    stop("Unknown formula for y_slope")

}
### <---------------------------------------------------------------------->
#' @rdname ecld.y_slope
ecld.y_slope_trunc <- function(object, t=1) {
	
    ecld.validate(object, sged.allowed=TRUE)
    one <- if(object@use.mpfr) ecd.mp1 else 1 # for consistent type

    lambda <- object@lambda * one
    s <- object@sigma * one
    beta <- object@beta
    mu <- object@mu

    if (lambda <= 2) return(ecd.mpnum(object,Inf))
    
    # SGED
    if (object@is.sged) {
        sp = s*(1+beta)
        slope <- sp * (2/(sp*t*lambda))^(lambda/(lambda-2))
        return(slope+object@mu)
	}

    # symmetric
    if (beta==0) {
        slope <- s * (2/(s*t*lambda))^(lambda/(lambda-2))
        return(slope+object@mu)
    }

    # general
    
    slope_eq <- function(x) {
        y <- ecld.solve(object, x)
        p <- 2*(x-mu)/s + beta*y
        q <- -lambda*s*(-y)^(lambda-1) - beta*(x-mu)
        p/q+t
    }

    lower <- 0.1*s + mu
    if (slope_eq(lower) > 0) {
        # no solution
        return(NaN * one)
    }
    
    # determine the lower/upper range
    upper <- 10 + mu
    repeat {
        if (slope_eq(upper/2) < 0) {
            lower <- (upper-mu)/2 + mu
        }
        if (slope_eq(upper) > 0) break
        upper <- (upper-mu)*10 + mu
        
    }
    
    if (object@use.mpfr) {
        nmax <- unirootR(slope_eq, lower=lower, upper=upper)
        return(nmax$root)
    } else {
        nmax <- uniroot(slope_eq, lower=lower, upper=upper)
        return(nmax$root)
    }

    stop("Unknown formula for y_slope_trunc")

}