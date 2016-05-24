#' Slope of \eqn{y(x)}
#' 
#' Slope of \eqn{y(x)}, that is, \eqn{dy/dx}. 
#'
#' @method y_slope ecd
#'
#' @param object an object of ecd class
#' @param x a numeric vector of \code{x} dimension
#'
#' @return a numeric vector of \eqn{dy/dx}
#'
#' @keywords elliptic-curve
#'
#' @author Stephen H. Lihn
#'
#' @export y_slope
#'
#' @examples
#' d <- ecd(0,1)
#' x <- seq(-20,20,by=0.01)
#' yp <- y_slope(d,x)

### <======================================================================>
"y_slope.ecd" <- function(object, x)
{
    N <- length(x)
    
    a <- object@alpha
    r <- object@gamma
    s <- object@sigma
    b <- object@beta
    mu <- object@mu
    L <- object@lambda
    
    f <- function(x) {
        xi <- (x-mu)/s
        y <- solve(object, x)
        yL <- L*(-y)^(L-1) # for L=3, this is just 3*y^2
        -(b*y+2*xi)/(yL+b*xi+r)/s
    }
    rs <- simplify2array(lapply(x, f))
    return(ecd.mpnum(object,rs))
}
### <---------------------------------------------------------------------->
#' @rdname y_slope.ecd
setGeneric("y_slope", function(object, x) standardGeneric("y_slope"))
#' @rdname y_slope.ecd
setMethod("y_slope", signature("ecd"), y_slope.ecd)
### <---------------------------------------------------------------------->
