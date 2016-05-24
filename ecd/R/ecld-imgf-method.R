#' Incomplete moment generating function (IMGF) of ecld
#'
#' The analytic solutions for IMGF of ecld, if available.
#' Note that, by default, risk neutrality is honored. However, you must note that
#' when fitting market data, this is usually not true.
#' SGED is supported.
#'
#' @param object an object of ecld class
#' @param k a numeric vector of log-strike
#' @param otype character, specifying option type:
#'              \code{c} (default) or \code{p}.
#' @param RN logical, use risk-neutral assumption for \code{mu_D}
#'
#' @return numeric, incomplete MGF
#'
#' @keywords mgf
#'
#' @author Stephen H-T. Lihn
#'
#' @export ecld.imgf
#' @export ecld.imgf_gamma
#' @export ecld.imgf_integrate
#'
#' @examples
#' ld <- ecld(sigma=0.01)
#' ecld.imgf(ld,0)
#'
### <======================================================================>
"ecld.imgf" <- function(object, k, otype="c", RN=TRUE)
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
    M1 <- if (RN) 1 else ecld.mgf(object)
    
    # SGED
    if (object@is.sged) {
        s2 <- ecld.ifelse(object, k<mu, s*(1-b), s*(1+b))
        ki <- (k-mu)/s # not adjust for sigma +/-
        
        # summation
        nmax <- floor(ecd.mp2f(ecld.mgf_trunc(object)))
        if (is.na(nmax)) {
            stop("NA found in mgf_trunc! Is sigma too large?")
        }
        if (nmax > 100) nmax <- 100 # cap the length of summation

        otype2 <- function(ki) {
            if (ki>=0 & otype=="p") return("c")
            if (ki<0 & otype=="c") return("p")
            return(otype)
        }

        sum_terms <- function(ki) {
            n <- seq(0, nmax) # use all terms, no skip
            f <- function(n) ecld.imnt(object, ki, order=n, otype=otype2(ki))/gamma(n+1)
            terms <- ecld.sapply(object, n, f)
            if (any(is.na(terms))) {
                stop("NA found in mgf_term. Consider using MPFR to improve precision!")
            }
            sum(terms)
        }
        
        calc_imgf <- function(ki) {
            G2 <- ecd.mpnum(object, exp(mu)*sum_terms(ki))
            if (otype==otype2(ki)) return(ecd.mpnum(object, G2))
            else return(ecd.mpnum(object, M1-G2))
        }
        return(ecld.sapply(object, ki, calc_imgf))
    }

    # normal
    if (lambda==1) {
        if (b != 0) {
            stop("lambda=1: beta must be zero")
        }
        if (RN) {
            m <- 1/2*ecd.erf(k/s-s/4)
            if (otype=="c") return(1/2-m)
            if (otype=="p") return(1/2+m)
        } else {
            mg <- ecld.mgf(object, t=1)
            ki <- (k-object@mu)/s
            m <- ecd.erf(ki-s/2)
            if (otype=="c") return(1/2*mg*(1-m))
            if (otype=="p") return(1/2*mg*(1+m))
        }
    }

    ki <- (k-mu)/s
	
    if (lambda==2) {
        sgn <- ifelse(ki<0, -1, 1) # sign of k
	    B0 <- ecld.laplace_B(b, 0)
	    Bs <- ecld.laplace_B(b, -sgn, s) # sigma extension
	    y <- -Bs*abs(ki)
	    m <- 1/2/Bs/B0 * exp(y) * exp(mu)
	    
	    Mc <- ifelse(ki>=0, m, M1-m)
	    Mp <- ifelse(ki<0,  m, M1-m)
	    if (otype=="c") return(ecd.mpnum(object, Mc))
	    if (otype=="p") return(ecd.mpnum(object, Mp))
	}
	
    # The remaining parametrization must integrate directly
    m <- ecld.imgf_integrate(object, k, otype=otype, RN=RN)
    return(m)
    
    stop("Unknown analytic formula for IMGF")

}
### <---------------------------------------------------------------------->
#' @rdname ecld.imgf
"ecld.imgf_gamma" <- function(object, k, otype="c", RN=TRUE)
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
    M1 <- if (RN) 1 else ecld.mgf(object)
    
    ki <- (k-mu)/s
    
    # SGED
    if (object@is.sged) {
        s <- ecld.ifelse(object, k<mu, s*(1-b), s*(1+b)) # override s
        ki <- (k-mu)/s # adjusted for sigma +/-
    } else if (lambda==2) {
        
        sgn_ki <- ifelse(ki<0, -1, 1) # sign of ki
        B <- ecld.laplace_B(b, -sgn_ki)
        B0 <- ecld.laplace_B(b, 0)
        ki <- (k-mu)/s*B # adjusted for sigma +/-
        order <- 100
        mnt2 <- function(n) {
            s^n * B^(-n-1)/B0 * ecld.gamma(n+1, abs(ki))
        }
        im <- mnt2(0)
        for (n in 1:order) {
            sgn <- ecd.ifelse(object, ki>=0, 1, (-1)^n)
            im <- im + sgn*mnt2(n)/gamma(n+1)
        }
        im <- im * exp(mu)/2
        
        Mc <- ifelse(ki>=0, im, M1-im)
        Mp <- ifelse(ki<0,  im, M1-im)
        if (otype=="c") return(ecd.mpnum(object, Mc))
        if (otype=="p") return(ecd.mpnum(object, Mp))
        stop("Unknown option")
        
        
    } else if (object@beta != 0) {
        stop("Beta must be zero")
    }
    
    order <- if (object@is.sged) ecld.mgf_trunc(object) else ecld.y_slope_trunc(object)
    if (order > 100) order <- 100
    
    k2 <- abs(ki)^(2/lambda)
    mnt <- function(n) s^n * ecld.gamma(lambda*(n+1)/2, k2)/gamma(lambda/2)
    
    im <- mnt(0)
    for (n in 1:floor(order)) {
        sgn <- ecd.ifelse(object, ki>=0, 1, (-1)^n)
        im <- im + sgn*mnt(n)/gamma(n+1)
    }
    im <- im * exp(mu)/2
    
    if (object@is.sged) {
        im <- im * s / object@sigma
    }
    
    Mc <- ifelse(ki>=0, im, M1-im)
    Mp <- ifelse(ki<0,  im, M1-im)
    if (otype=="c") return(ecd.mpnum(object, Mc))
    if (otype=="p") return(ecd.mpnum(object, Mp))
    
    stop("Unknown option")
}
### <---------------------------------------------------------------------->
#' @rdname ecld.imgf
"ecld.imgf_integrate" <- function(object, k, otype="c", RN=TRUE)
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
    M1 <- if (RN) 1 else ecld.mgf(object)
    ki <- (k-mu)/s
    
    
    if (length(k) > 1) {
        # TODO This is okay, but could be better!
        f <- function(k) ecld.imgf_integrate(object, k, otype=otype, RN=RN)
        M <-  simplify2array(parallel::mclapply(k, f))
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
            e_y_n <- function(x) exp(-x^(2/lambda) - sn*x)
            M <- ecd.integrate(d0, e_y_n, -ki, Inf)
        } else {
            xt <- ecld.y_slope_trunc(object)/s2
            if (xt > .ecd.mpfr.N.sigma) xt <- .ecd.mpfr.N.sigma
            e_y_p <- function(x) exp(-x^(2/lambda) + sp*x)
            M <- ecd.integrate(d0, e_y_p, ki, xt)
        }
        if (M$message != "OK") {
            stop("Failed to integrate SGED IMGF from unit distribution")
        }
        m <- M$value*s2/ecld.const(object)*exp(mu)
        
        Mc <- ifelse(ki>=0, m, M1-m)
        Mp <- ifelse(ki<0,  m, M1-m)
        if (otype=="c") return(ecd.mpnum(object, Mc))
        if (otype=="p") return(ecd.mpnum(object, Mp))
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
    e_y <- function(xi) exp(ecld.solve(ld0,xi) + s*xi)
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
    
    Mc <- ifelse(ki>=0, m, M1-m)
    Mp <- ifelse(ki<0,  m, M1-m)
    if (otype=="c") return(ecd.mpnum(object, Mc))
    if (otype=="p") return(ecd.mpnum(object, Mp))

    stop("Unknown option")
    
}
### <---------------------------------------------------------------------->
