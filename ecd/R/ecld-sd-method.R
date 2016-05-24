#' Compute statistics analytically for an ecld object
#' 
#' Compute statistics for mean, var, skewness, kurtosis,
#' from the known analytical result. SGED is supported.
#'
#' @param object an object of ecld class
#'
#' @return numeric or mpfr
#'
#' @keywords statistics
#'
#' @author Stephen H-T. Lihn
#'
#' @export ecld.sd
#' @export ecld.var
#' @export ecld.mean
#' @export ecld.skewness
#' @export ecld.kurt
#' @export ecld.kurtosis
#'
#' @examples
#' ld <- ecld(3)
#' ecld.sd(ld)
#' ecld.var(ld)
#' ecld.mean(ld)
#' ecld.skewness(ld)
#' ecld.kurt(ld)
#'
### <======================================================================>
"ecld.sd" <- function(object)
{
    sqrt(ecld.var(object))
}
### <---------------------------------------------------------------------->
#' @rdname ecld.sd
"ecld.var" <- function(object)
{
    ecld.validate(object, sged.allowed=TRUE)
    s <- object@sigma
    l <- object@lambda

    # SGED
    if (object@is.sged) {
        # use moment
        v <- ecld.moment(object,2)-ecld.moment(object,1)^2
        return(v)
    }
    
    # symmetric
    if (object@beta==0) {
        x <- gamma(l*3/2)
        y <- gamma(l/2)
        return(s^2*x/y)
    }
    
    # use moment
    v <- ecld.moment(object,2)-ecld.moment(object,1)^2
    return(v)
    
    stop("Unknown analytic formula for var")
}
### <---------------------------------------------------------------------->
#' @rdname ecld.sd
"ecld.mean" <- function(object)
{
    ecld.validate(object, sged.allowed=TRUE)

    l <- object@lambda

    # SGED
    if (object@is.sged) {
        G <- gamma(l) * gamma(l/2)
        v <- 2*object@sigma *object@beta *G
        return(v)
    }

    if (object@beta==0) {
        return(object@mu)
    }
    stop("Unknown analytic formula for mean")
}
### <---------------------------------------------------------------------->
#' @rdname ecld.sd
"ecld.skewness" <- function(object)
{
    ecld.validate(object)
    s <- object@sigma
    l <- object@lambda
    
    if (object@beta==0) {
        return(0)
    }
    stop("Unknown analytic formula for skewness")
}
### <---------------------------------------------------------------------->
#' @rdname ecld.sd
"ecld.kurtosis" <- function(object)
{
    ecld.validate(object)
    s <- object@sigma
    l <- object@lambda
    
    if (object@beta==0) {
        x <- gamma(l/2) * gamma(l*5/2)
        y <- gamma(l*3/2)^2
        return(x/y)
    }
    stop("Unknown analytic formula for kurtosis")
}
### <---------------------------------------------------------------------->
#' @rdname ecld.sd
"ecld.kurt" <- function(object) ecld.kurtosis(object)
### <---------------------------------------------------------------------->



