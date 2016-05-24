#' Integrate a function with PDF of the distribution
#' 
#' Integrate a function with PDF of the distribution. 
#' The integration is seperated into three segments to ensure convergence.
#'
#' @method integrate_pdf ecd
#'
#' @param object An object of ecd class
#' @param f An R function taking a numeric first argument and 
#'          returning a numeric vector of the same length. 
#'          Returning a non-finite element will generate an error.
#' @param lower	Numeric, the lower limit of integration. Can be infinite.
#' @param upper Numeric, the upper limit of integration. Can be infinite.
#' @param ... Addtional arguments for \code{f}.
#' @param show.warning logical, display warning messages.
#' @param verbose logical, display timing information, for debugging purpose.
#'
#' @return A list of class "\code{integrate}".
#'
#' @keywords integrate pdf
#'
#' @author Stephen H. Lihn
#'
#' @export integrate_pdf
#'
#' @importFrom Rmpfr integrateR
#' @importFrom stats integrate
#'
#' @examples
#' d <- ecd()
#' integrate_pdf(d, function(x){x^2}, -Inf, Inf)
### <======================================================================>
"integrate_pdf.ecd" <- function(object, f, lower, upper, ...,
                                show.warning=TRUE, verbose=FALSE)
{
    # debug on bisection
    if (verbose) print(paste(Sys.time(), "integrate_pdf",
                        "lower=", ecd.mp2f(lower),
                        "upper=", ecd.mp2f(upper)
                        ))
    
    if (length(lower)>1) {
        stop("parameter 'lower' should be length-one!")
    }
    if (length(upper)>1) {
        stop("parameter 'upper' should be length-one!")
    }
    if (is.na(lower)) {
        stop(paste("parameter 'lower' cannot be NA! upper=", ecd.mp2f(upper)))
    }
    if (is.na(upper)) {
        stop(paste("parameter 'upper' cannot be NA! lower=", ecd.mp2f(lower)))
    }
    
    # ------------------------------------------------- #
    # make sure lower < upper
    if (lower == upper) return(.ecd.intg_zero(object))
    if (lower > upper) {
        p <- integrate_pdf(object, f, upper, lower, ..., show.warning=show.warning)
        p$value <- -p$value
        return(p)
    }
    
    # ------------------------------------------------- #
    # now lower < upper is always true
 
    ipdf <- function(lower2, upper2, ...) {
        if (lower2 == upper2) return(.ecd.intg_zero(object))

        p <- integrate_pdf(object, f, lower2, upper2, ...,
                           show.warning=show.warning, verbose=verbose)
        
        if (show.warning & p$message != "OK") {
            warning(paste("WARN: integrate_pdf", p$message,
                          "lower=", ecd.mp2f(lower2), "upper=", ecd.mp2f(upper2),
                          "for ecd:", ecd.toString(object)))
        }
        
        p
    }

    # split integral WRT c1 and c2
    sigma <- object@sigma
    c1 <- object@const_left_x
    c2 <- object@const_right_x

    if (lower < c1-sigma & upper > c1) {
        p1 <- ipdf(lower, c1, ...)
        p2 <- ipdf(c1, upper, ...)
        p2$abs.error <- p1$abs.error + p2$abs.error
        p2$value <- p1$value + p2$value
        return(p2)
    }
    if (upper > c2+sigma & lower < c2) {
        p1 <- ipdf(lower, c2, ...)
        p2 <- ipdf(c2, upper, ...)
        p2$abs.error <- p1$abs.error + p2$abs.error
        p2$value <- p1$value + p2$value
        return(p2)
    }
    
    # mpfr, to further break it down by 4*sd beyond c1 and c2
    if (object@use.mpfr) {
        sd <- object@R^(1/3)*object@sigma
        d1 <- object@mu - 4*sd
        d2 <- object@mu + 4*sd

        if (lower < c1 & lower < d1-sigma & upper > d1) {
            p1 <- ipdf(lower, d1, ...)
            p2 <- ipdf(d1, upper, ...)
            p2$abs.error <- p1$abs.error + p2$abs.error
            p2$value <- p1$value + p2$value
            return(p2)
        }
        if (upper > c2 & upper > d2+sigma & lower < d2) {
            p1 <- ipdf(lower, d2, ...)
            p2 <- ipdf(d2, upper, ...)
            p2$abs.error <- p1$abs.error + p2$abs.error
            p2$value <- p1$value + p2$value
            return(p2)
        }
    }

	# odd cases
    if (lower == -Inf & upper == -Inf) return(.ecd.intg_zero(object))
    if (lower == Inf  & upper == Inf)  return(.ecd.intg_zero(object))

	# otherwise, do the normal thing

    # expected value of f
    E_of_f <- function(x)
	{
    	ecd.pdf(object, x)*f(x, ...)
	}

	# debug, when integration blows up, especially for imgf
    if (verbose) print(paste(Sys.time(), "call ecd.integrate:",
                      "lower=", ecd.mp2f(lower),
                      "upper=", ecd.mp2f(upper)
                      ))

    ecd.integrate(object, E_of_f, lower, upper, show.warning=show.warning)
}
### <---------------------------------------------------------------------->
#' @rdname integrate_pdf.ecd
setGeneric("integrate_pdf", function(object, f, lower, upper, ...) standardGeneric("integrate_pdf"))
#' @rdname integrate_pdf.ecd
setMethod("integrate_pdf", signature("ecd"), integrate_pdf.ecd)
### <---------------------------------------------------------------------->
.ecd.intg_zero <- function(object) {
	if (object@use.mpfr) {
        Rmpfr::integrateR(function(x){x},0,0)
    } else {
        integrate(function(x){x},0,0)
    }
}
### <---------------------------------------------------------------------->
