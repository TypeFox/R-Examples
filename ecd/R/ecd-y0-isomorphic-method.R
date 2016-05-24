#' The analytic solution of \eqn{y(0)} via isomorphic mapping.
#' 
#' This utility can be called two ways: 
#' (a) specify R and theta; (b) provide the ecd object.
#' But not at the same time.
#'
#' @param theta numeric vector, the polar coordinate
#' @param R numeric vector, the polar coordinate
#' @param object optionally, a single ecd object
#'
#' @return the value of y(0)
#'
#' @keywords solve
#'
#' @export
#'
#' @examples
#' t <- 45/180*pi
#' ecd.y0_isomorphic(t)
#'
### <======================================================================>
"ecd.y0_isomorphic" <- function(theta=NaN, R=1, object=NULL)
{
    if (!is.null(object)) {
        if (!is.na(theta)) {
            stop("Can not specify theta and object at the same time!")
        }
        theta <- object@theta
        R <- object@R
    }
    
    # handle vector
    if (length(theta)>1 | length(R)>1) {
        if (length(theta)==1) {
            theta <- rep(0,length(R)) + theta
        }
        if (length(R)==1) {
            R <- rep(0,length(theta)) + R
        }
        rs <- theta*0
        for(i in 1:length(theta)) {
            rs[i] <- ecd.y0_isomorphic(theta[i], R[i])
        }
        return(rs)
    }
    
    # handle single computation
    
    # cusp
    if (abs(theta/pi-7/4) <= .Machine$double.eps) {
        alpha <- R*cos(theta)
        return(-(alpha/2)^(1/3))
    }
    # j=0, positive
    if (theta == 0) {
        return(R^(1/3))
    }
    
    fs <- -2^(2/3) *abs(sin(theta))^(1/3) *R^(1/3)
    if (theta <= pi & theta > 0) {
        T <- log(tan(theta/2))/3
        return(fs*sinh(T))
    }
    if (theta > pi & theta < 7/4*pi) {
        if(theta > 5/4*pi) theta <- theta*(1+0i)
        # theta needs to be real between degree (180, 225]
        # theta needs to be complex between degree (225, 315)
        T2 <- acosh(1/tan(theta))/3
        y0 <- fs*cosh(T2)
        
        if (is.na(y0)) return(Re(y0))
        if (abs(Im(y0)) < 1e-6) return(Re(y0))
        
        # this should only occur past the critical line
        warning(paste("Unexpected imaginary result:", y0))
        return(y0)
        
    }
    stop(paste("theta out of range:", theta))
}
### <---------------------------------------------------------------------->
