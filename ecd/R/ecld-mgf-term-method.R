#' The term structure of ecld symmetric MGF
#'
#' \code{ecld.mgf_term} and \code{ecld.mgf_diterm} are the term and derivative
#' of the term by order (n) in the summation of MGF.
#' \code{ecld.mgf_trunc} uses \code{ecld.mgf_diterm} to locate the truncation
#' of MGF terms.
#' \code{ecld.mgf_trunc_max_sigma} locates the maximum sigma that keeps MGF finite for each lambda.
#' SGED is supported.
#'
#' @param object an object of ecd class
#' @param order numeric. Order of the term (moment)
#' @param t numeric, for MGF
#'
#' @return numeric
#'
#' @keywords moment
#'
#' @author Stephen H-T. Lihn
#'
#' @export ecld.mgf_term
#' @export ecld.mgf_diterm
#' @export ecld.mgf_trunc
#' @export ecld.mgf_trunc_max_sigma
#'
#' @importFrom stats uniroot
#'
#' @examples
#' ld <- ecld(3, sigma=0.01*ecd.mp1)
#' ecld.mgf_trunc(ld)
#'
### <======================================================================>
"ecld.mgf_term" <- function(object, order, t=1)
{
    ecld.validate(object, sged.allowed=TRUE)
    one <- if(object@use.mpfr) ecd.mp1 else 1 # for gamma function

    lambda <- object@lambda * one
    s <- object@sigma * one
    mu <- object@mu * one
    n <- order * one
    
    # SGED
    if (object@is.sged) {
        # use moment
        v <- ecld.moment(object,n) *t^n /gamma(n+1)
        return(v * exp(t*mu))
    }
    
    if (object@beta==0) {
        x <- gamma(lambda*(n+1)/2)
        y <- gamma(lambda/2)*gamma(n+1)
        return((s*t)^n * x/y * exp(t*mu))
    }
    stop("Unknown analytic formula for MGF term")

}
### <---------------------------------------------------------------------->
#' @rdname ecld.mgf_term
"ecld.mgf_diterm" <- function(object, order, t=1)
{
    ecld.validate(object, sged.allowed=TRUE)
    one <- if(object@use.mpfr) ecd.mp1 else 1 # for gamma function
    
    lambda <- object@lambda * one
    s <- object@sigma * one
    b <- object@beta * one
    n <- order * one

    # SGED
    if (object@is.sged) {
        b <- object@beta * ecd.mp1 # need MPFR for very large n
        dlogFdn <- if (n*log(1+abs(b)) > 20) log(1+abs(b)) else ecld.diterm_sged_dlogF_dn(n, b)
        dGdn <- log(s*t) + digamma(lambda*(n+1)/2) *lambda/2 - digamma(n+1)
        di <- dGdn + dlogFdn
        if(! object@use.mpfr) di <- ecd.mp2f(di)
        return(di)
    }

    if (object@beta==0) {
        dGdn <- log(s*t) + digamma(lambda*(n+1)/2) *lambda/2 - digamma(n+1)
        return(dGdn)
    }
    stop("Unknown analytic formula for MGF diterm")
    
}
# the following is a utility, it is not exported, but it is tested
# dlogF_dn -> log(1+abs(b)), as n->Inf
"ecld.diterm_sged_dlogF_dn" <- function(n, b) {
    P <- function(n, b) (1+b)^n
    F <- function(n, b) (P(n,b) + P(n,-b))/2
    H <- function(n, b) (log(1+b)*P(n,b) + log(1-b)*P(n,-b))/2
    dlogFdn <- H(n+1,b)/F(n+1,b)
    return(dlogFdn)
}
### <---------------------------------------------------------------------->
#' @rdname ecld.mgf_term
"ecld.mgf_trunc" <- function(object, t=1) {
    ecld.validate(object, sged.allowed=TRUE)
    one <- if(object@use.mpfr) ecd.mp1 else 1 # for gamma function
    
    lambda <- object@lambda * one
    s <- object@sigma
    n <- order

    # lambda <= 2, there is no need for truncation
    if (lambda <= 2) {
        return(Inf)
    }
    diterm <- function(x) {
        di <- ecld.mgf_diterm(object, x, t=t)
        if(is.na(di)) {
            stop(paste("NA is found in diterm: x=", x, " for ",
                "lambda=", ecd.mp2f(lambda),
                "sigma=", ecd.mp2f(object@sigma),
                "beta=", ecd.mp2f(object@beta),
                "is.sged=", object@is.sged))
        }
        return(di)
    }
    
    lower <- 0.1*s
    if (diterm(lower) > 0) {
        # no solution
        return(NaN * one)
    }
    
    # determine the lower/upper range
    upper <- 10
    repeat {
        if (diterm(upper/2) < 0) {
            lower <- upper/2
        }
        if (diterm(upper) > 0) break
        upper <- upper*10
        
    }
    
    # root finding
    if (object@use.mpfr) {
        nmax <- unirootR(diterm, lower=lower, upper=upper)
        nmax$root
    } else {
        nmax <- uniroot(diterm, lower=lower, upper=upper)
        nmax$root
    }
}
### <---------------------------------------------------------------------->
#' @rdname ecld.mgf_term
"ecld.mgf_trunc_max_sigma" <- function(object, order=1) {
    ecld.validate(object)
    one <- if(object@use.mpfr) ecd.mp1 else 1 # for gamma function
    
    lambda <- object@lambda * one
    n <- order * one
    exp(digamma(n+1) - lambda/2*digamma(lambda/2*(n+1)))
}
### <---------------------------------------------------------------------->


