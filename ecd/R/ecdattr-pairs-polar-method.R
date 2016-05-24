#' Create a list of basic \code{ecdattr} objects in polar coordinate
#' 
#' The list is created by the Cartesian product between \code{R} and \code{theta}.
#' This contains the data points of a circular area defined by \code{R,theta}.
#' If \code{cusp} is 1, data points are on the critical line specified by \code{R}.
#'
#' @param R,theta numeric vectors
#' @param cusp numeric, representing type of cusp. Only 0 (default) and 1 are allowed.
#' @param use.mpfr logical, whether to use mpfr for ecd object, default is \code{FALSE}.
#'
#' @return a list of basic \code{ecdattr} objects.
#'
#' @keywords ecdattr
#'
#' @export
#'
### <======================================================================>
ecdattr.pairs_polar <- function(R, theta, cusp=0, use.mpfr=FALSE) {
    lst <- list()
    cnt <- 1
    if (cusp == 1) {
        if (! is.nan(theta)) stop("theta must be NaN when cusp=1")
        theta <- 315/180*pi
        alpha <- R*cos(theta)
        for (a in alpha) {
            lst[[as.character(cnt)]] <- ecdattr(a, NaN, cusp, use.mpfr)
            cnt <- cnt+1
        }
        return(lst)
    } else if (cusp == 0) {
        for (r in R) {
            for (t in theta) {
                lst[[as.character(cnt)]] <- ecdattr(r*cos(t), r*sin(t), cusp, use.mpfr)
                cnt <- cnt+1
            }
        }
    } else {
        stop("Cusp must be 0 or 1")
    }
    lst
}
### <---------------------------------------------------------------------->
