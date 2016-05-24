#' Incomplete moment (imnt) of ecld
#'
#' The analytic solutions for imnt of ecld, if available.
#' Note that, by default, risk neutrality is honored.
#' \code{ecld.imnt_sum} provides an alternative method to calculate IMGF.
#'
#' @param object an object of ecld class
#' @param ki numeric vector of normalized log-strike, \code{(k-mu)/sigma}
#' @param order numeric. Order of the moment to be computed.
#'              For \code{ecld.imnt_sum}, this is the maximum order to be truncated.
#'              For small sigma at lambda=3, this can be simply 2.
#'              If \code{Inf}, the slope truncation procedure will be used
#'              to determine the maximum order. However, due to the numeric
#'              limit of \code{pgamma}, it is capped at 100.
#' @param otype character, specifying option type:
#'              \code{c} (default) or \code{p}.
#'
#' @return numeric vector
#'
#' @keywords mgf
#'
#' @author Stephen H-T. Lihn
#'
#' @export ecld.imnt
#' @export ecld.imnt_integrate
#' @export ecld.imnt_sum
#'
#' @examples
#' ld <- ecld(sigma=0.01*ecd.mp1)
#' ki <- seq(-0.1, 0.1, by=0.01)
#' ecld.imnt(ld,ki, 1)
### <======================================================================>
"ecld.imnt" <- function(object, ki, order, otype="c")
{
    if (!(otype %in% c("c","p"))) {
        stop(paste("Unknown option type:", otype))
    }

    ecld.validate(object, sged.allowed=TRUE)
    one <- if(object@use.mpfr) ecd.mp1 else 1 # for gamma function
    
    lambda <- object@lambda * one
    s <- object@sigma * one
    b <- object@beta
    n <- order
    
    # SGED
    if (object@is.sged) {
        n_sgn <- ifelse(ki>=0, 1, (-1)^n)
        s2 <- ecld.ifelse(object, ki<0, s*(1-b), s*(1+b))
        kii <- ki*s/s2 # adjust for sigma +/-
        k2 <- abs(kii)^(2/lambda)
        G2 <- n_sgn * ecld.gamma(lambda*(n+1)/2, k2)/gamma(lambda/2)/2 * s2^n
        # s2^n: multiplier from unit to actual
        G2 <- ecd.mpnum(object, G2*s2/s)
        
        M1 <- ecld.moment(object, n)
        Mc <- ifelse(ki>=0, G2, M1-G2)
        Mp <- ifelse(ki<0,  G2, M1-G2)
        
        if (otype=="c") return(ecd.mpnum(object, Mc))
        if (otype=="p") return(ecd.mpnum(object, Mp))
    }

    # symmetric
    if (b == 0) {
        n_sgn <- ifelse(ki>=0, 1, (-1)^n)
        k2 <- abs(ki)^(2/lambda)
        G2 <- n_sgn * ecld.gamma(lambda*(n+1)/2, k2)/gamma(lambda/2)/2 * s^n
        # s^n: multiplier from unit to actual
        G2 <- ecd.mpnum(object, G2)
        
        M1 <- ecld.moment(object, n)
        Mc <- ifelse(ki>=0, G2, M1-G2)
        Mp <- ifelse(ki<0,  G2, M1-G2)
        
        if (otype=="c") return(ecd.mpnum(object, Mc))
        if (otype=="p") return(ecd.mpnum(object, Mp))
    }
	
    if (lambda==2) {
        sgn <- ifelse(ki<0, 1, -1)
	    B0 <- ecld.laplace_B(b, 0)
	    B <- ecld.laplace_B(b, sgn)
        n_sgn <- ifelse(ki<0, (-1)^n, 1)
	    m <- 1/2/B0/B^(n+1) * ecld.gamma(n+1, abs(ki)*B) * n_sgn * s^n

	    M1 <- ecld.moment(object, n)
	    Mc <- ecd.ifelse(object, ki>=0, m, M1-m)
	    Mp <- ecd.ifelse(object, ki<0,  m, M1-m)
	    if (otype=="c") return(Mc)
	    if (otype=="p") return(Mp)
	}
	
    # The remaining parametrization must integrate directly
    m <- ecld.imnt_integrate(object, ki, order, otype=otype)
    return(m)
    
    stop("Unknown analytic formula for imnt")

}
### <---------------------------------------------------------------------->
#' @rdname ecld.imnt
"ecld.imnt_integrate" <- function(object, ki, order, otype="c")
{
    if (!(otype %in% c("c","p"))) {
        stop(paste("Unknown option type:", otype))
    }
    
    ecld.validate(object)
    one <- if(object@use.mpfr) ecd.mp1 else 1 # for gamma function
    
    lambda <- object@lambda * one
    s <- object@sigma * one
    b <- object@beta
    n <- order
    
    if (length(ki) > 1) {
        # TODO This is okay, but could be better!
        f <- function(ki) ecld.imnt_integrate(object, ki, order, otype=otype)
        M <-  simplify2array(parallel::mclapply(ki, f))
        #M <-  NaN*ki 
        #for (i in 1: length(ki)) {
        #    M[i] <- f(ki[i])
        #}
        return(ecd.mpnum(object, M))
    }

    # now handle length-one ki
    M1 <- ecld.moment(object, n)

    # TODO
    #if (ki==Inf) return(ecd.mpnum(object, 1))
    #if (ki==-Inf) return(ecd.mpnum(object, 0))
    
    # MPFR is channelled through sigma=1
    # since we are using unit distribution, either way should be fine
    ld0 <- ecld(lambda=ecd.mp2f(lambda), beta=ecd.mp2f(b), sigma=one)
    d0 <- ecd(lambda=ecd.mp2f(lambda), beta=ecd.mp2f(b), sigma=one, bare.bone=TRUE)
    M <- NULL
    e_y <- function(xi) exp(ecld.solve(ld0,xi)) * abs(xi)^n
    if (ki < 0) {
        M <- ecd.integrate(d0, e_y, -Inf, ki)
    } else {
        M <- ecd.integrate(d0, e_y, ki, Inf)
    }
    if (M$message != "OK") {
        stop("Failed to integrate moment from unit distribution")
    }
    n_sgn <- ifelse(ki<0, (-1)^n, 1)
    m <- M$value/ecld.const(ld0) * n_sgn * s^n
    
    Mc <- ecd.ifelse(object, ki>=0, m, M1-m)
    Mp <- ecd.ifelse(object, ki<0,  m, M1-m)
    if (otype=="c") return(Mc)
    if (otype=="p") return(Mp)
    
    stop("Unknown option")
    
}
### <---------------------------------------------------------------------->
#' @rdname ecld.imnt
"ecld.imnt_sum" <- function(object, ki, order, otype="c")
{
    if (!(otype %in% c("c","p"))) {
        stop(paste("Unknown option type:", otype))
    }
    
    ecld.validate(object)
    one <- if(object@use.mpfr) ecd.mp1 else 1 # for gamma function

    if (order == Inf) {
        order <- ecld.y_slope_trunc(object)
        if (is.na(order)) {
            stop(paste("max order of sum truncation (y_slope_trunc) is NaN:",
                 "lambda=", ecd.mp2f(object@lambda),
                 "beta=", ecd.mp2f(object@beta),
                 "sigma=", ecd.mp2f(object@sigma)
                ))
        }
    }
    stopifnot(!is.na(order) & order >= 0)
    if (order > 100) order <- 100
    
    im <- 0*ki # init
    for (n in 0:floor(order)) {
        im <- im + ecld.imnt(object, ki, n, otype=otype)/gamma(n+1)
    }
    # this result does NOT include exp(mu)
    return(ecd.mpnum(object, im))
}
### <---------------------------------------------------------------------->
