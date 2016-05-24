#' CDF and CCDF of ecld
#'
#' The analytic solutions for CDF and CCDF of ecld, if available.
#' \code{ecld.cdf_gamma} is a sub-module with the CDF expressed as
#' incomplete gamma function.
#' SGED is supported only in \code{ecld.cdf} and \code{ecld.ccdf}.
#'
#' @param object an object of ecld class
#' @param x a numeric vector of \code{x}
#'
#' @return The CDF or CCDF vector
#'
#' @keywords cdf
#'
#' @author Stephen H. Lihn
#'
#' @export ecld.cdf
#' @export ecld.ccdf
#' @export ecld.cdf_gamma
#' @export ecld.cdf_integrate
#'
#' @examples
#' ld <- ecld(sigma=0.01*ecd.mp1)
#' x <- seq(-0.1, 0.1, by=0.01)
#' ecld.cdf(ld,x)
### <======================================================================>
"ecld.cdf" <- function(object, x)
{
    ecld.validate(object, sged.allowed=TRUE)
    one <- if(object@use.mpfr) ecd.mp1 else 1 # for gamma function
    
    lambda <- object@lambda * one
    s <- object@sigma * one
    b <- object@beta
    mu <- object@mu
	xi <- (x-mu)/s
	
    # SGED
    if (object@is.sged) {
        s2 <- ecld.ifelse(object, x<mu, s*(1-b), s*(1+b))
        xi <- (x-mu)/s2
        g <- ecld.gamma(lambda/2, abs(xi)^(2/lambda))
        c <- g /gamma(lambda/2) *s2/s/2
        return(ecld.ifelse(object, x<mu, c, 1-c))
    }
    
    # normal
    if (lambda==1) {
        if (b != 0) {
            stop("lambda=1: beta must be zero")
        }
        return(1/2*(ecd.erf(xi)+1))
    }
    
    if (lambda==2) {
        sgn <- ifelse(xi<0, -1, 1) # sign of x
    	B0 <- ecld.laplace_B(b, 0)
        B <- ecld.laplace_B(b, -sgn)
        cd <- 1/2/B/B0 * exp(-B*abs(xi))
        return(ecd.mpnum(object, ifelse(xi<0, cd, 1-cd)))
    }

	if (lambda==3 & b==0) {
	    p <- 1/sqrt(pi)*abs(xi)^(1/3)*exp(-abs(xi)^(2/3))
	    q <- 1/2*(1-ecd.erf(abs(xi)^(1/3)))
	    return(ecd.mpnum(object, ifelse(xi<0, p+q, 1-p-q)))
	}
    
    # for all other sym distributions
    if (b==0) {
        return(ecld.cdf_gamma(object, x))
    }

    # for asym dist with fractional lambda, must integrate to get CDF
	if (lambda > 2) {
        return(ecld.cdf_integrate(object, x))
	}
	
    stop("Unknown analytic formula for CDF")

}
### <---------------------------------------------------------------------->
#' @rdname ecld.cdf
"ecld.ccdf" <- function(object, x)
{
    1-ecld.cdf(object, x)
}
### <---------------------------------------------------------------------->
#' @rdname ecld.cdf
"ecld.cdf_integrate" <- function(object, x)
{
    ecld.validate(object, sged.allowed=TRUE)
    one <- if(object@use.mpfr) ecd.mp1 else 1 # for gamma function
    
    lambda <- object@lambda * one
    s <- object@sigma * one
    b <- object@beta
    mu <- object@mu
	xi <- (x-mu)/s
	   
    if (length(x) > 1) {
        # TODO This is okay, but could be better!
        f <- function(x) ecld.cdf_integrate(object,x)
        M <-  simplify2array(parallel::mclapply(x, f))
        return(ecd.mpnum(object, M))
    }
    
    # SGED
    if (object@is.sged) {
        sp <- s*(1+b)
        sn <- s*(1-b)
        s2 <- ecld.ifelse(object, x<mu, sn, sp)
        xi <- (x-mu)/s2
        
        d0 <- ecd(lambda=ecd.mp2f(lambda), beta=ecd.mp2f(b), sigma=one, bare.bone=TRUE)
        e_y <- function(x) exp(-x^(2/lambda))
        M <- ecd.integrate(d0, e_y, abs(xi), Inf)
        if (M$message != "OK") {
            stop("Failed to integrate SGED IMGF from unit distribution")
        }
        C <- ecld.const(object)/s2
        
        if (xi < 0) {
            return(ecd.mpnum(object, M$value/C))
        } else {
            return(ecd.mpnum(object, 1-M$value/C))
        }
    }
    
    if (xi==Inf) return(ecd.mpnum(object, 1))
    if (xi==-Inf) return(ecd.mpnum(object, 0))
    
    # MPFR is channelled through sigma=1
    # since we are using unit distribution, either way should be fine
    ld0 <- ecld(lambda=ecd.mp2f(lambda), beta=ecd.mp2f(b), sigma=one)
    d0 <- ecd(lambda=ecd.mp2f(lambda), beta=ecd.mp2f(b), sigma=one, bare.bone=TRUE)
    e_y0 <- function(xi) exp(ecld.solve(ld0,xi))
    M <- NULL
    if (xi < 0) {
        M <- ecd.integrate(d0, e_y0, -Inf, xi)
    } else {
        M <- ecd.integrate(d0, e_y0, xi, Inf)
    }
    if (M$message != "OK") {
        stop("Failed to integrate moment from unit distribution")
    }
    C <- ecld.const(ld0)
    
    if (xi < 0) {
        return(ecd.mpnum(object, M$value/C))
    } else {
        return(ecd.mpnum(object, 1-M$value/C))
    }
    
}
### <---------------------------------------------------------------------->
#' @rdname ecld.cdf
"ecld.cdf_gamma" <- function(object, x)
{
    ecld.validate(object)
    one <- if(object@use.mpfr) ecd.mp1 else 1 # for gamma function

    lambda <- object@lambda * one
    s <- object@sigma * one
    xi <- (x-object@mu)/s
    
    if (object@beta != 0) {
        stop("This is CDF for symmetric distribution. Beta must be zero")
    }

    x2 <- abs(xi)^(2/lambda)
    G2 <- ecld.gamma(lambda/2, x2) / gamma(lambda/2) /2
    return(ecd.mpnum(object, ifelse(xi <= 0, G2, 1-G2)))

}
### <---------------------------------------------------------------------->
