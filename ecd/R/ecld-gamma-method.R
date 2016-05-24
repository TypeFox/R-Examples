#' Incomplete gamma function and asymptotic expansion
#'
#' \code{ecld.gamma} is the wrapper for incomplete gamma function
#' \eqn{\Gamma(s,x)}. It is mainly to wrap around \code{pgamma}.
#' And \code{ecld.gamma_hgeo} is the asymptotic expansion of \eqn{\Gamma(s,x)}
#' using hypergeometric series, \eqn{e^{-x} x^{s-1} {}_2 F_0 (1,1-s;;-1/x)}.
#' It is mainly used in for star OGF \eqn{L^{*}(k;\lambda)}.
#' \code{ecld.gamma_2F0} is simply \eqn{{}_2 F_0 (1,1-s;;-1/x)}, which is used
#' in the star OGF expansion.
#'
#' @param s numeric
#' @param x numeric
#' @param order numeric, the order of the power series
#'
#' @return numeric
#'
#' @keywords gamma
#'
#' @author Stephen H-T. Lihn
#'
#' @export ecld.gamma
#' @export ecld.gamma_hgeo
#' @export ecld.gamma_2F0
#'
#' @importFrom stats pgamma
#'
### <======================================================================>
"ecld.gamma" <- function(s, x=0) {
    # TODO pgamma can't handle MPFR
    # TODO s seems to have a limit of about 100-200
    pgamma(ecd.mp2f(x), ecd.mp2f(s), lower.tail = FALSE)*gamma(s)
}
### <---------------------------------------------------------------------->
#' @rdname ecld.gamma
"ecld.gamma_hgeo" <- function(s, x, order) {
    p <- exp(-x) * x^(s-1)
    q <- ecld.gamma_2F0(s, x, order)
    return(p*q)
}
### <---------------------------------------------------------------------->
#' @rdname ecld.gamma
"ecld.gamma_2F0" <- function(s, x, order) {
    if (order == 0) return(1)
    stopifnot(order > 0)
    
    q <- 1
    xi <- 1
    si <- 1
    for (i in 1:order) {
        si <- si * (s-i)
        xi <- xi * x
        q <- q + si/xi
    }
    return(q)
}
### <---------------------------------------------------------------------->

