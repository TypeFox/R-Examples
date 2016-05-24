#' Star OGF of ecld
#'
#' The star OGF of ecld is the limiting OGF for small sigma.
#' It only depends on the normalized k and lambda.
#' Its dependency on sigma and mu is removed.
#' SGED is not supported yet.
#'
#' @param object an object of ecld class
#' @param ki a numeric vector of log-strike
#' @param order numeric, order of the hypergeometric series to be computed
#'
#' @return The state price of option in star OGF terms.
#'
#' @keywords ogf
#'
#' @author Stephen H-T. Lihn
#'
#' @export ecld.ogf_star
#' @export ecld.ogf_star_hgeo
#' @export ecld.ogf_star_exp
#' @export ecld.ogf_star_gamma_star
#'
#' @examples
#' ld <- ecld(sigma=0.001*ecd.mp1)
#' ki <- seq(1, 5, by=1)
#' ecld.ogf_star(ld, ki)
### <======================================================================>
"ecld.ogf_star" <- function(object, ki)
{
    ecld.validate(object)
    one <- if(object@use.mpfr) ecd.mp1 else 1 # for gamma function
    
    lambda <- object@lambda * one
    xi <- abs(ki)^(2/lambda)
    
    q <- 1/2/gamma(lambda/2)
    p <- ecld.gamma(lambda, xi) - abs(ki)*ecld.gamma(lambda/2, xi)
    return(p*q)
}
### <---------------------------------------------------------------------->
#' @rdname ecld.ogf_star
"ecld.ogf_star_hgeo" <- function(object, ki, order=4)
{
    ecld.validate(object)
    one <- if(object@use.mpfr) ecd.mp1 else 1 # for gamma function
    
    lambda <- object@lambda * one
    xi <- abs(ki)^(2/lambda)

    q <- exp(-xi)*abs(ki)^(2-2/lambda)/2/gamma(lambda/2)
    p <- ecld.gamma_2F0(lambda, xi, order) - ecld.gamma_2F0(lambda/2, xi, order)
    return(p*q)
}
### <---------------------------------------------------------------------->
#' @rdname ecld.ogf_star
"ecld.ogf_star_exp" <- function(object, ki, order=3)
{
    ecld.validate(object)
    one <- if(object@use.mpfr) ecd.mp1 else 1 # for gamma function
    
    l <- object@lambda * one
    xi <- abs(ki)^(2/l) 

    q <- exp(-abs(ki)^(2/l))*abs(ki)^(2-4/l) *l/4/gamma(l/2)
    if (order == 0) return(q)
    p <- 1
    if (l==1 & order >= 1) {
        dfac <- 1
        for (n in 1:order) {
            dfac <- dfac*(2*n+1)
            p <- p + dfac/(-2*one*xi)^n
        } 
        return (p*q)
    }
    if (order >= 1) p <- p + 3/2*(l-2)/xi
    if (order >= 2) p <- p + 1/4*(7*l^2 -36*l +44)/xi^2
    if (order >= 3) p <- p + 1/8*(15*l^3 -140*l^2 +420*l -400)/xi^3
    if (order >= 4) {
        stop("order >= 4 is not support at the moment")
    }
    return(p*q)
}
### <---------------------------------------------------------------------->
#' @rdname ecld.ogf_star
"ecld.ogf_star_gamma_star" <- function(object, ki, order=6)
{
    ecld.validate(object)
    one <- if(object@use.mpfr) ecd.mp1 else 1 # for gamma function
    
    lambda <- object@lambda * one
    
    k2l <- abs(ki)^(2/lambda)
    q <- gamma(lambda)/gamma(lambda/2)/2 - abs(ki)/2
    
    p <- 0
    for (n in 0:order) {
        f1 = 1/gamma(lambda/2+n+1)
        f2 = gamma(lambda)/gamma(lambda/2)/gamma(lambda+n+1)
        p <- p + k2l^n*(f1-f2)
    }
    return(q + ki^2/2*exp(-k2l)*p)
}
### <---------------------------------------------------------------------->
