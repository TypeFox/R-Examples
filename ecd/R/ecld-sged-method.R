#' The integral solutions of SGED
#' 
#' These integrals are mainly used as validation to analytic solutions.
#' If you must use them, be mindful of their slower speeds.
#'
#' @param object an sged object of ecld class
#' @param x a numeric vector of \code{x}
#' @param order numeric, order of the moment to be computed
#' @param k a numeric vector of log-strike
#' @param t numeric, for MGF and IMGF
#' @param otype character, specifying option type:
#'              \code{c} (default) or \code{p}.
#'
#' @return numeric 
#'
#' @keywords sged
#'
#' @author Stephen H-T. Lihn
#'
#' @export ecld.sged_const
#' @export ecld.sged_cdf
#' @export ecld.sged_moment
#' @export ecld.sged_mgf
#' @export ecld.sged_imgf
#' @export ecld.sged_ogf
#'
#' @examples
#' ld <- ecld(3)
#' ecld.const(ld)

### <======================================================================>
"ecld.sged_const" <- function(object)
{
    ecld.validate(object, sged.only=TRUE)
    one <- if(object@use.mpfr) ecd.mp1 else 1 # for gamma function

    lambda <- object@lambda
    s <- object@sigma
    b <- object@beta
    mu <- object@mu
    
    d <- ecd(lambda=lambda, sigma=s, mu=mu, bare.bone=TRUE)
    
    e_y <- function(x) exp(ecld.solve(object, x))
    i1 <- ecd.integrate(d, e_y, mu, Inf)
    i2 <- ecd.integrate(d, e_y, -Inf, mu)
    C <- i1$value + i2$value
    return(C)
}
### <---------------------------------------------------------------------->
#' @rdname ecld.sged_const
"ecld.sged_cdf" <- function(object, x)
{
    ecld.validate(object, sged.only=TRUE)
    
    d <- ecd(lambda=object@lambda, sigma=object@sigma, bare.bone=TRUE)
    e_y <- function(x) exp(ecld.solve(object, x))
    C <- ecld.const(object)
    
    cdf <- function(x) {
        if (x < object@mu) {
            c2 <- ecd.integrate(d, e_y, -Inf, x)
            return(c2$value/C)
        } else {
            c1 <- ecd.integrate(d, e_y, x, Inf)
            return(1-c1$value/C)
        }
    }
    return(ecld.sapply(object, x, cdf))
}
### <---------------------------------------------------------------------->
#' @rdname ecld.sged_const
"ecld.sged_moment" <- function(object, order)
{
    # always ignore mu, always central moment
    ecld.validate(object, sged.only=TRUE)
    
    if (length(order)>1) {
        return(ecld.sapply(object, order, function(n) ecld.sged_moment(object,n) ))
    }
    
    d <- ecd(lambda=object@lambda, sigma=object@sigma, bare.bone=TRUE)
    e_y <- function(x) exp(ecld.solve(object, x)) * (x-object@mu)^order
    C <- ecld.const(object)
    c1 <- ecd.integrate(d, e_y, object@mu, Inf)
    c2 <- ecd.integrate(d, e_y, -Inf, object@mu)
    return((c1$value + c2$value)/C)
}
### <---------------------------------------------------------------------->
#' @rdname ecld.sged_const
"ecld.sged_mgf" <- function(object, t=1)
{
    ecld.validate(object, sged.only=TRUE)
    
    sigma <- object@sigma
    d <- ecd(lambda=object@lambda, sigma=sigma, bare.bone=TRUE)
    e_y <- function(x) exp(ecld.solve(object, x) + t*x)
    C <- ecld.const(object)
    
    xmax <- ecld.y_slope_trunc(object)
    if (is.na(xmax)) {
        stop("Failed to locate y_slope truncation point")
    }
    xmax2 <- .ecd.mpfr.N.sigma * sigma + object@mu
    if (xmax > xmax2) xmax <- xmax2
    
    c1 <- ecd.integrate(d, e_y, object@mu, xmax)
    c2 <- ecd.integrate(d, e_y, -Inf, object@mu)
    return(c1$value/C + c2$value/C)
}
### <---------------------------------------------------------------------->
#' @rdname ecld.sged_const
"ecld.sged_imgf" <- function(object, k, t=1, otype="c")
{
    # This function integrates directly, so it uses object@mu
    if (length(k)>1) {
        f <- function(k) ecld.sged_imgf(object, k, t=t, otype=otype)
        return(ecld.sapply(object, k, f))
    }
    
    if (!(otype %in% c("c","p"))) {
        stop(paste("Unknown option type:", otype))
    }

    ecld.validate(object, sged.only=TRUE)
    
    d <- ecd(lambda=object@lambda, sigma=object@sigma, bare.bone=TRUE)
    e_y <- function(x) exp(ecld.solve(object, x) + t*x)
    C <- ecld.const(object)
    
    if (otype=="c") {
        if (k >= object@mu) {
            xmax <- ecld.y_slope_trunc(object)
            if (is.na(xmax)) {
                stop("Failed to locate y_slope truncation point")
            }
            xmax2 <- .ecd.mpfr.N.sigma * object@sigma + object@mu
            if (xmax > xmax2) xmax <- xmax2
            c1 <- ecd.integrate(d, e_y, k, xmax)
            return(c1$value/C)
        } else {
            Mp <- ecld.sged_imgf(object, k, otype="p")
            M1 <- ecld.sged_mgf(object)
            return(M1-Mp)
        }
    }
    if (otype=="p") {
        if (k < object@mu) {
            c2 <- ecd.integrate(d, e_y, -Inf, k)
            return(c2$value/C)
        } else {
            Mc <- ecld.sged_imgf(object, k, otype="c")
            M1 <- ecld.sged_mgf(object)
            return(M1-Mc)
        }
    }
    
    stop(paste("Unknown option type:", otype))
}
### <---------------------------------------------------------------------->
#' @rdname ecld.sged_const
"ecld.sged_ogf" <- function(object, k, otype="c")
{
    # This function integrates directly, so it uses object@mu
    if (length(k)>1) {
        f <- function(k) ecld.sged_ogf(object, k, otype=otype)
        return(ecld.mclapply(object, k, f))
    }
    
    if (!(otype %in% c("c","p"))) {
        stop(paste("Unknown option type:", otype))
    }
    
    ecld.validate(object, sged.only=TRUE)
    
    # can't use sged_cdf, it causes too much error 
    if (otype=="c") {
        Mc <- ecld.imgf(object, k, otype="c")
        ccdf <- 1-ecld.cdf(object, k)
        return(Mc-exp(k)*ccdf)
    }
    if (otype=="p") {
        Mp <- ecld.imgf(object, k, otype="p")
        cdf <- ecld.cdf(object, k)
        return(-Mp+exp(k)*cdf)
    }
    
    stop(paste("Unknown option type:", otype))
}
### <---------------------------------------------------------------------->

