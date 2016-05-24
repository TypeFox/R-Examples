#' Ellipticity of ecd object
#' 
#' Ellipticity of ecd object, defined as half of the distance
#' between the two elliptic points.
#'
#' @method ellipticity ecd
#'
#' @param object An object of ecd class
#' @param tol Numeric, the tolerance of precision during subdivision.
#'            Default: \code{1e-4} of stdev.
#'
#' @return a list with 3 major numbers: xe1= negative x_e, xe2= positive x_e, avg= ellipticity
#'
#' @keywords ellipticity
#'
#' @export ellipticity
#'
#' @examples
#' d <- ecd(0,1)
#' ellipticity(d)
### <======================================================================>
"ellipticity.ecd" <- function(object, tol=1e-4)
{
    # handle singularity for cusp
    if (object@cusp > 0) {
        if (object@alpha == 0 & object@gamma == 0) {
            return(list(xe1=0, xe2=0, avg=0, dx=c(0,0), iter=c(1,1)))
        }
    }
    
    if (length(object@stats)==0) {
        stop("stats is not computed in ecd object")
    }

    stdev <- ecd.mp2f(object@stats$stdev)
    m1 <- ecd.mp2f(object@stats$mean)

    N <- 3 # starting with N stdev range
    M <- 1.25 # subdivide M * dx on each side in an iteration

    f <- function(init.x.from, init.x.to, direction) {
        iter <- 1
        xe <- .ecd.find_elliptic_point(object, init.x.from, init.x.to, direction, 100)
        xe$iter <- iter
        repeat {
        	dx <- xe$dx
            xe <- .ecd.find_elliptic_point(object, xe$x-M*dx, xe$x+M*dx, direction, 10)
            iter <- iter+1
            if (is.na(dx) | is.na(stdev) | stdev==0) break
            if (dx/stdev < tol | iter >= 12) break
        }
        xe$iter <- iter
        xe
    }
    
    xe1 <- f(m1-N*stdev, m1, -1)
    xe2 <- f(m1, m1+N*stdev, 1)
    avg <- (xe2$x-xe1$x)/2
    list(xe1=xe1$x, xe2=xe2$x, avg=avg, dx=c(xe1$dx, xe2$dx), iter=c(xe1$iter, xe2$iter))
}
### <---------------------------------------------------------------------->
#' @rdname ellipticity.ecd
setGeneric("ellipticity", function(object, tol=1e-5) standardGeneric("ellipticity"))
#' @rdname ellipticity.ecd
setMethod("ellipticity", signature("ecd"), ellipticity.ecd)
### <---------------------------------------------------------------------->
# helper function
.ecd.find_elliptic_point <- function(object, x.from, x.to, direction, N) {
    if (is.na(x.from) | is.na(x.to)) {
        return(list(x=NaN, dx=NaN))
    }
    
    dx <- (x.to-x.from)/N
    x <- seq(x.from, x.to, by=dx)
    yp <- ecd.mp2f(y_slope(object, x))
    if (direction > 0) {
        x2 <- min(x[yp==min(yp)])
        return(list(x=x2, dx=dx))
    } else {
        x1 <- max(x[yp==max(yp)])
        return(list(x=x1, dx=dx))
    }
}