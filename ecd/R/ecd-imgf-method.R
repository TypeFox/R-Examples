#' Incomplete MGF of ecd
#'
#' Incomplete moment generating function (IMGF) of ecd, integration of \eqn{e^z P(z)} for z from \code{x} to \code{Inf}.
#' \code{ecd.mu_D} is simply a wrapper around MGF.
#'
#' @param object an object of ecd class
#' @param x a numeric vector of \code{x}, default to \code{-Inf}
#' @param t a numeric value for MGF, default to \code{1}
#' @param minus1 logical, subtracting one from \eqn{e^{tx}}
#' @param unit.sigma logical, transforming to unit sigma to achieve greater stability.
#'                    Due to the instability of quadpack for \code{ecd.integrate_pdf}, default to \code{TRUE}.
#'                    But constructing a new ecd object has significant overhead, be aware of it in performance sensitive program.
#' @param n.sigma length-one numeric, specifying the max number of sigma to check for truncation.
#' @param verbose logical, display timing information, for debugging purpose.
#'
#' @return The IMGF
#'
#' @keywords cdf option
#'
#' @author Stephen H. Lihn
#'
#' @export ecd.imgf
#' @export ecd.mu_D
#'
#' @examples
#' d <- ecd(0, 0, sigma=0.01)
#' x <- seq(0, 1, by=0.1)
#' ecd.imgf(d, x)
### <======================================================================>
"ecd.imgf" <- function(object, x=-Inf, t=1, minus1=FALSE, unit.sigma=FALSE, n.sigma=.ecd.mpfr.N.sigma, verbose=FALSE)
{
    if (length(object@stats)==0) {
        stop("stats is not computed in ecd object")
    }

    s <- object@sigma
    m <- object@stats$mean

    if (verbose) print(paste(Sys.time(), "ecd.imgf: start",
                            "x=", ecd.mp2f(x), "t=", ecd.mp2f(t)
                     ))

    if (unit.sigma & object@sigma != 1) {
	    object2 <- ecd(alpha= object@alpha, 
	        		   gamma= object@gamma,
	  				   sigma= 1,
	  				   beta=  object@beta,
	  				   mu=    object@mu/s)
	    I <- ecd.imgf(object2, x=x/s, t=s, unit.sigma=FALSE, verbose=verbose)
	    return(I)
	}

    # figure out the truncation
	fn_yt <- function(x) solve(object,x) + t*x
    xs <- m+s*seq(1, n.sigma, by=5)
    yt <- fn_yt(xs)
    max_x <- max(xs[yt==min(yt)])
    stopifnot(is.numericMpfr(max_x))
    if (verbose) print(paste(Sys.time(), "ecd.imgf: use max_x=", ecd.mp2f(max_x) ))
    
    # execute the integral
    x <- ecd.mpnum(object, ifelse(x <= max_x, x, max_x))
    c <- ifelse(minus1, 1, 0)
    fn_exp <- function(z) exp(t*z)-c
	imgf <- ecd.ccdf(object, x, to.x = max_x, f = fn_exp, piece.wise=TRUE, verbose=verbose)

    if (verbose) print(paste(Sys.time(), "ecd.imgf: done. imgf=", ecd.mp2f(imgf) ))
    imgf
}
### <---------------------------------------------------------------------->
#' @rdname ecd.imgf
"ecd.mu_D" <- function(object)
{
    return(-log(ecd.imgf(object)))
}

### <---------------------------------------------------------------------->
