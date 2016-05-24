#' Wrapper to integrate numeric and mpfr
#'
#' The wrapper handles chooses to to use \code{integrate} for numeric;
#' or to use \code{integrateR} for mpfr. Since the later doesn't allow
#' infinity, there is a special handling to replace infinity with a large
#' multiple of \code{sigma}.
#'
#' @param object An object of ecd class. This object can be bare-boned.
#' @param f An R function taking a numeric first argument and
#'          returning a numeric vector of the same length.
#'          Returning a non-finite element will generate an error.
#' @param lower	Numeric, the lower limit of integration. Can be infinite.
#' @param upper Numeric, the upper limit of integration. Can be infinite.
#' @param abs.tol numeric, the suggested absolute tolerance.
#' @param ... Addtional arguments for \code{f}.
#' @param mpfr.qagi logical, to use quadpack qagi transformation for infinity.
#' @param show.warning logical, to suppress warnings or not.
#'
#' @return The \code{integrate} object
#'
#' @keywords utility
#'
#' @export
#'
#' @importFrom stats integrate
#'
#' @author Stephen H. Lihn
#'
#' @importFrom Rmpfr integrateR
#'
### <======================================================================>
ecd.integrate <- function(object, f, lower, upper, ...,
                          abs.tol=.Machine$double.eps^0.25,
                          mpfr.qagi=TRUE, show.warning=TRUE) {

    # double precision
	if (!object@use.mpfr) {
        return(integrate(f, lower, upper, ..., abs.tol=abs.tol, stop.on.error=FALSE))
    }

    # ------------------------------------------------ #
    # mpfr mode
    has.inf <- abs(lower)+abs(upper) == Inf
    
    # make sure upper and lower are in the right order
    if (has.inf & lower > upper) {
        return(ecd.integrate(object, function(x) -f(x), upper, lower, ...,
                             abs.tol=abs.tol, mpfr.qagi=mpfr.qagi,
                             show.warning=show.warning))
    }
    
    # qagi
    if (has.inf & mpfr.qagi) {
        R <- ecd.mp2f(object@R)
        degree <- floor(ecd.mp2f(object@theta/pi*180))
        tzero <- ecd.mpfr(.Machine$double.eps)
        intg3 <- function(fn, a, b) ecd.integrate(object, fn, a, b, ...,
                                                  abs.tol=abs.tol, mpfr.qagi=TRUE,
                                                  show.warning=show.warning)
        if (lower < 0 & upper > 0) {
            p1 <- intg3(f, lower, 0)
            p2 <- intg3(f, 0, upper)
            return(.ecd.intg_merge(p1, p2))
        }
        if ((lower == -Inf & upper <= 0) | (upper == Inf & lower >= 0)) {
            return(ecd.mpfr_qagi(object, f, lower, upper, ...,
                                 abs.tol=abs.tol, show.warning=show.warning))
        }
        stop("ecd.integrate: mpfr.qagi: Unhandled condition?")
        
    } else if (has.inf & !mpfr.qagi) {
        # truncation method, very slow, but not dependent on my newly developed qagi method
        N <- .ecd.mpfr.N.sigma # number of sigma as replacement for +/- Inf
        if (lower == -Inf) {
            lower <- -N*object@sigma
        }
        if (upper == Inf) {
            upper <- N*object@sigma
        }
    }

    # finally, execute it
    intg <- function() {
        abs.tol <- ecd.adj_abs_tol(f, lower, upper, abs.tol)
        f2 <- function(x) f(x,...)
        Rmpfr::integrateR(f2, lower, upper, abs.tol=abs.tol)
    }
    if (show.warning) intg() else suppressWarnings(intg())

}
### <---------------------------------------------------------------------->
.ecd.intg_merge <- function(p1, p2) {
    p1$value <- p1$value + p2$value
    p1$abs.error <- p1$abs.error + p2$abs.error
    p1$subdivisions <- p1$subdivisions + p2$subdivisions
    
    if (p2$message != "OK") {
        if(p1$message == "OK") {
            p1$message <- p2$message
        } else {
            p1$message <- paste(p1$message, " <and> ", p2$message)
        }
    }
    return(p1)
}
### <---------------------------------------------------------------------->

