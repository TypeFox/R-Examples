#' Create a list of basic \code{ecdattr} objects
#' 
#' The list is created by the Cartesian product between \code{alpha} and \code{gamma}.
#' This contains the data points of a rectangular area defined by \code{alpha,gamma}.
#' If \code{cusp} is 1, data points are on the critical line specified by \code{alpha}.
#'
#' @param alpha,gamma numeric vectors
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
ecdattr.pairs <- function(alpha, gamma, cusp=0, use.mpfr=FALSE) {
    lst <- list()
    cnt <- 1
    if (cusp == 1) {
        if (! is.nan(gamma)) stop("gamma must be NaN when cusp=1")
        for (a in alpha) {
            lst[[as.character(cnt)]] <- ecdattr(a, NaN, cusp, use.mpfr)
            cnt <- cnt+1
        }
        return(lst)
    } else if (cusp == 0) {
        for (a in alpha) {
            for (r in gamma) {
                lst[[as.character(cnt)]] <- ecdattr(a, r, cusp, use.mpfr)
                cnt <- cnt+1
            }
        }
    } else {
        stop("Cusp must be 0 or 1")
    }
    lst
}
### <---------------------------------------------------------------------->
