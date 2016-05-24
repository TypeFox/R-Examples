#' Option generating function (OGF) of ecld
#'
#' The analytic solutions for OGF of ecld, if available.
#' Note that, by default, risk neutrality is honored. However, you must note that
#' when fitting market data, this is usually not true.
#'
#' @param object an object of ecld class
#' @param k a numeric vector of log-strike
#' @param order numeric, order of the moment to be computed
#' @param otype character, specifying option type:
#'              \code{c} (default) or \code{p}.
#' @param RN logical, use risk-neutral assumption for \code{mu_D}
#'
#' @return The state price of option
#'
#' @keywords ogf
#'
#' @author Stephen H-T. Lihn
#'
#' @export ecld.ogf
#' @export ecld.ogf_integrate
#' @export ecld.ogf_gamma
#' @export ecld.ogf_imnt_sum
#' @export ecld.ogf_log_slope
#'
#' @examples
#' ld <- ecld(sigma=0.01*ecd.mp1)
#' k <- seq(-0.1, 0.1, by=0.05)
#' ecld.ogf(ld,k)
### <======================================================================>
"ecld.ogf" <- function(object, k, otype="c", RN=TRUE)
{
    if (!(otype %in% c("c","p"))) {
        stop(paste("Unknown option type:", otype))
    }

    ecld.validate(object, sged.allowed=TRUE)
    one <- if(object@use.mpfr) ecd.mp1 else 1 # for gamma function
    
    lambda <- object@lambda * one
    s <- object@sigma * one
    b <- object@beta
	
    # SGED, use ogf_gamma
    if (object@is.sged) {
        m <- ecld.ogf_gamma(object, k, otype=otype, RN=RN)
        return(m)
    }
    
    # normal
    if (lambda==1) {
        if (b != 0) {
            stop("lambda=1: beta must be zero")
        }
        if (! RN) {
            stop("lambda=1: RN must be true")
        }

        M1k <- 1-exp(k)
        p <- -1/2*ecd.erf(k/s-s/4)
        q <- exp(k)/2*ecd.erf(k/s+s/4)
        Lc <- p + q + M1k/2
        if (otype=="c") return(Lc)
        if (otype=="p") return(Lc - M1k)
    }

    mu <- if (RN) ecld.mu_D(object) else object@mu
    M1 <- if (RN) 1 else exp(object@mu - ecld.mu_D(object))
    M1k <- M1-exp(k)
    ki <- (k-mu)/s
    
    if (lambda==2) {
        sgn <- ecd.ifelse(object, ki<0, -1, 1) # sign of k
        B <- function(b, sgn=0, s=0) ecld.laplace_B(b, sgn, s)
	    B0 <- B(b, 0)
	    p <- 1/B(b, -sgn, s) - 1/B(b, -sgn)
	    q <- exp(-B(b, -sgn, s)*abs(ki))
	    m <- 1/2/B0 * exp(mu) * p * q
        # m is one-sided Lc or -Lp
        
	    Lc  <- ecd.ifelse(object, ki>=0, m, M1k-m)
	    NLp <- ecd.ifelse(object, ki<0,  m, M1k-m)
	    if (otype=="c") return(ecd.mpnum(object, Lc))
	    if (otype=="p") return(ecd.mpnum(object, -NLp))
    }
    
    # symmetric case
    if (b==0) {
        m <- ecld.ogf_gamma(object, k, otype=otype, RN=RN)
        return(m)
    }

    # The remaining parametrization must integrate directly
    m <- ecld.ogf_integrate(object, k, otype=otype, RN=RN)
    return(m)

    stop("Unknown analytic formula for OGF")

}
### <---------------------------------------------------------------------->
#' @rdname ecld.ogf
"ecld.ogf_integrate" <- function(object, k, otype="c", RN=TRUE)
{
    if (!(otype %in% c("c","p"))) {
        stop(paste("Unknown option type:", otype))
    }
    
    ecld.validate(object, sged.allowed=TRUE)
    one <- if(object@use.mpfr) ecd.mp1 else 1 # for gamma function
    
    lambda <- object@lambda * one
    s <- object@sigma * one
    b <- object@beta

    mu <- if (RN) ecld.mu_D(object) else object@mu
    M1 <- if (RN) 1 else exp(object@mu - ecld.mu_D(object))
    M1k <- M1-exp(k)
    ki <- (k-mu)/s
    
    
    if (length(k) > 1) {
        # TODO This is okay, but could be better!
        f <- function(k) ecld.ogf_integrate(object, k, otype=otype, RN=RN)
        M <- simplify2array(parallel::mclapply(k, f))
        return(ecd.mpnum(object, M))
    }
    
    # SGED
    if (object@is.sged) {
        sp <- s*(1+b)
        sn <- s*(1-b)
        s2 <- ecld.ifelse(object, k<mu, sn, sp)
        ki <- (k-mu)/s2
        
        d0 <- ecd(lambda=ecd.mp2f(lambda), beta=ecd.mp2f(b), sigma=one, bare.bone=TRUE)
        M <- NULL
        if (ki < 0) {
            e_y_n <- function(x) exp(-x^(2/lambda)) *(exp(-sn*x) - exp(sn*ki))
            M <- ecd.integrate(d0, e_y_n, -ki, Inf)
        } else {
            xt <- ecld.y_slope_trunc(object)/s2
            if (xt > .ecd.mpfr.N.sigma) xt <- .ecd.mpfr.N.sigma
            e_y_p <- function(x) exp(-x^(2/lambda)) * (exp(sp*x) - exp(sp*ki))
            M <- ecd.integrate(d0, e_y_p, ki, xt)
        }
        if (M$message != "OK") {
            stop("Failed to integrate SGED IMGF from unit distribution")
        }
        m <- M$value * s2/ecld.const(object)*exp(mu)
        
        Lc <- ecd.ifelse(object, ki>=0, m, M1k-m)
        NLp <- ecd.ifelse(object, ki<0, m, M1k-m)
        if (otype=="c") return(ecd.mpnum(object, Lc))
        if (otype=="p") return(ecd.mpnum(object, -NLp))
        stop("Unknown otype")
    }

    # TODO
    #if (ki==Inf) return(ecd.mpnum(object, 1))
    #if (ki==-Inf) return(ecd.mpnum(object, 0))
    
    # MPFR is channelled through sigma=1
    # since we are using unit distribution, either way should be fine
    ld0 <- ecld(lambda=ecd.mp2f(lambda), beta=ecd.mp2f(b), sigma=one)
    d0 <- ecd(lambda=ecd.mp2f(lambda), beta=ecd.mp2f(b), sigma=one, bare.bone=TRUE)
    M <- NULL
    e_y <- function(xi) exp(ecld.solve(ld0,xi)) * (exp(s*xi) - exp(s*ki))
    if (ki < 0) {
        M <- ecd.integrate(d0, e_y, -Inf, ki)
    } else {
        xt <- ecld.y_slope_trunc(object)/s
        if (xt > .ecd.mpfr.N.sigma) xt <- .ecd.mpfr.N.sigma
        M <- ecd.integrate(d0, e_y, ki, xt)
    }
    if (M$message != "OK") {
        stop("Failed to integrate IMGF from unit distribution")
    }
    m <- M$value/ecld.const(ld0)*exp(mu)
    
    Lc  <- ecd.ifelse(object, ki>=0, m, M1k-m)
    NLp <- ecd.ifelse(object, ki<0,  m, M1k-m)
    if (otype=="c") return(ecd.mpnum(object, Lc))
    if (otype=="p") return(ecd.mpnum(object, -NLp))
    
    stop("Unknown option")
    
}
### <---------------------------------------------------------------------->
#' @rdname ecld.ogf
"ecld.ogf_gamma" <- function(object, k, otype="c", RN=TRUE)
{
    if (!(otype %in% c("c","p"))) {
        stop(paste("Unknown option type:", otype))
    }

    ecld.validate(object, sged.allowed=TRUE)
    one <- if(object@use.mpfr) ecd.mp1 else 1 # for gamma function

    lambda <- object@lambda * one
    s <- object@sigma * one
    b <- object@beta
    mu <- if (RN) ecld.mu_D(object) else object@mu
    M1 <- if (RN) 1 else exp(object@mu - ecld.mu_D(object))
    M1k <- M1-exp(k)
    ki <- (k-mu)/s

    # SGED
    if (object@is.sged) {
        s <- ecld.ifelse(object, k<mu, s*(1-b), s*(1+b)) # override s
        ki <- (k-mu)/s # adjusted for sigma +/-
    } else if (object@beta != 0) {
        stop("Beta must be zero")
    }

    order <- if (object@is.sged) ecld.mgf_trunc(object) else ecld.y_slope_trunc(object)
    if (order > 100) order <- 100

    k2 <- abs(ki)^(2/lambda)
    mnt <- function(n) s^n * ecld.gamma(lambda*(n+1)/2, k2)/gamma(lambda/2)

    im <- mnt(0)*(1-exp(s*ki))
    for (n in 1:floor(order)) {
        sgn <- ecd.ifelse(object, ki>=0, 1, (-1)^n)
        im <- im + sgn*mnt(n)/gamma(n+1)
    }
    im <- im * exp(mu)/2

    if (object@is.sged) {
        im <- im * s / object@sigma
    }
    
    Lc  <- ecd.ifelse(object, ki>=0, im, M1k-im)
    NLp <- ecd.ifelse(object, ki<0,  im, M1k-im)
    if (otype=="c") return(ecd.mpnum(object, Lc))
    if (otype=="p") return(ecd.mpnum(object, -NLp))
    
    stop("Unknown option")
}
### <---------------------------------------------------------------------->
#' @rdname ecld.ogf
"ecld.ogf_imnt_sum" <- function(object, k, order, otype="c", RN=TRUE)
{
    if (!(otype %in% c("c","p"))) {
        stop(paste("Unknown option type:", otype))
    }
    
    ecld.validate(object)
    one <- if(object@use.mpfr) ecd.mp1 else 1 # for gamma function

    mu <- if (RN) ecld.mu_D(object) else object@mu
    ki <- (k-mu)/object@sigma
    
    imgf <- ecld.imnt_sum(object, ki, order, otype=otype) * exp(mu)
    c <- exp(k) * ecld.imnt(object, ki, 0, otype=otype)
    L <- imgf-c
    if (otype=="c") return(ecd.mpnum(object, L))
    if (otype=="p") return(ecd.mpnum(object, -L))
    
    stop("Unknown option")
}
### <---------------------------------------------------------------------->
#' @rdname ecld.ogf
"ecld.ogf_log_slope" <- function(object, k, otype="c", RN=TRUE)
{
    if (!(otype %in% c("c","p"))) {
        stop(paste("Unknown option type:", otype))
    }
    
    ecld.validate(object)
    one <- if(object@use.mpfr) ecd.mp1 else 1 # for gamma function
    
    # construct RN object for CDF/CCDF
    object2 <- object
    if (RN) {
        mu <- ecld.mu_D(object)
        object2 <- ecld(lambda=object@lambda, sigma=object@sigma,
                        beta=object@beta, mu=mu)
    }
    cd <- if (otype=="c") -ecld.ccdf(object2, k) else ecld.cdf(object2, k)
    
    sd <- ecld.sd(object)
    p <- sd * exp(k) * cd
    q <- ecld.ogf(object, k, otype=otype, RN=RN)
    return(p/q)
}
### <---------------------------------------------------------------------->

