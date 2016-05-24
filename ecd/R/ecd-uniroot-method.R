#' Uniroot wrapper
#'
#' This function wraps ordinary uniroot and unirootR (from Rmpfr)
#' to the same interface.
#'
#' @param f the function for which the root is sought.
#' @param lower,upper the lower and upper end points of
#'                    the interval to be searched.
#' @param use.mpfr logical, to use MPFR (default), or else uniroot in stats.
#' @param tol the desired accuracy (convergence tolerance).
#' @param maxiter the maximum number of iterations.
#'
#' @return uniroot result
#'
#' @keywords utility
#'
#' @author Stephen H. Lihn
#'
#' @export
#'
#' @importFrom Rmpfr unirootR
#' @importFrom stats uniroot
#'
### <======================================================================>
"ecd.uniroot" <- function(f, lower, upper, use.mpfr = FALSE,
                          tol = .Machine$double.eps^0.25, maxiter = 1000)
{
    if (use.mpfr) {
        rt <- Rmpfr::unirootR(f, lower=lower, upper=upper,
                              tol=tol, maxiter=maxiter)
        return(rt)
    }
    else {
        rt <- uniroot(f, lower=ecd.mp2f(lower), upper=ecd.mp2f(upper),
                      tol=tol, maxiter=maxiter)
        return(rt)
    }
}
### <---------------------------------------------------------------------->
