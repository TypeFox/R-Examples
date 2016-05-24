#' Utility to integrate mpfr with infinity via qagi
#'
#' This utility supplements \code{Rmpfr::integrateR} with the quadpack qagi method to handle
#' integration involving infinity. Qagi is a transformation of \eqn{x/sigma=(1-t)/t} for positive x, 
#' and \eqn{x/sigma=(t-1)/t} for negative x. \eqn{t=0} is represented by \code{.Machine$double.eps}.
#' This utility requires (a) \code{lower} or \code{upper} is \code{+/-Inf}; 
#' (b) \code{lower} and \code{upper} are of the same sign.
#'
#' @param object an object of ecd class
#' @param f an R function taking a numeric first argument and
#'          returning a numeric vector of the same length.
#'          Returning a non-finite element will generate an error.
#' @param lower	numeric, the lower limit of integration. Can be infinite.
#' @param upper numeric, the upper limit of integration. Can be infinite.
#' @param abs.tol numeric, the suggested absolute tolerance.
#' @param ... addtional arguments for \code{f}.
#' @param show.warning logical, to suppress warnings or not.
#'
#' @return The \code{integrate} object
#'
#' @keywords utility
#'
#' @author Stephen H. Lihn
#'
#' @export
#'
#' @importFrom Rmpfr integrateR
#'
### <======================================================================>
ecd.mpfr_qagi <- function(object, f, lower, upper, ..., 
                          abs.tol=.Machine$double.eps^0.25, 
                          show.warning=TRUE) {

	if (!object@use.mpfr) {
        stop("Object doesn't use mpfr")
    }
    
    has.inf <- abs(lower)+abs(upper) == Inf
	same.sign <- sign(lower)*sign(upper)
    tzero <- ecd.mpfr(.Machine$double.eps)

    if (!has.inf) {
        stop("One of lower and upper must be infinity!")
    }
    if (same.sign<0) {
    	stop("Lower and upper must be of the same sign")
    }

    s <- object@sigma
    # ------------------------------------------------
	if (lower == -Inf & upper <= 0) {
		intg_minf <- function() {
			f_minf <- function(t) f((t-1)/t*s, ...)/t^2*s
			t <- 1/(1-upper/s)
			abs.tol <- ecd.adj_abs_tol(f_minf, tzero, t, abs.tol)
			Rmpfr::integrateR(f_minf, tzero, t, abs.tol=abs.tol)
		}
		return(if (show.warning) intg_minf() else suppressWarnings(intg_minf()))
	}        
    # ------------------------------------------------        
	if (upper == Inf & lower >= 0) {
		intg_inf <- function() {
			f_inf <- function(t) f((1-t)/t*s, ...)/t^2*s
			t <- 1/(1+lower/s)
			abs.tol <- ecd.adj_abs_tol(f_inf, tzero, t, abs.tol)
			Rmpfr::integrateR(f_inf, tzero, t, abs.tol=abs.tol)
		}
		return(if (show.warning) intg_inf() else suppressWarnings(intg_inf()))
	}
	stop("Unhandled condition?")
}
### <---------------------------------------------------------------------->
ecd.adj_abs_tol <- function(f, t1, t2, abs.tol, length.out=13, verbose=FALSE) {
    t <- t1 + (t2-t1)*seq(0, 1, length.out=length.out)
    dt <- (t2-t1)/length.out
    intg <- sum(abs(ecd.mpfr(sapply(t, f))))*dt
    abs.tol3 <- .Machine$double.eps^0.25*intg
    
    if (is.na(abs.tol3)) return(abs.tol)
    if (abs.tol3 > abs.tol) {
        abs.tol2 <- abs.tol
        abs.tol <- abs.tol3
        if (verbose) {
            print(paste("abs.tol is adjusted from", ecd.mp2f(abs.tol2), 
                        "to", ecd.mp2f(abs.tol),
                        "for intg=", ecd.mp2f(intg)
            ))
        }
    } else if (verbose) {
        print(paste("abs.tol stays at", ecd.mp2f(abs.tol), 
                    "for intg=", ecd.mp2f(intg)
             ))
    }
    ecd.mp2f(abs.tol)
}
