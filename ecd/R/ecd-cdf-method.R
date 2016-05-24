#' CDF of ecd
#'
#' CDF of ecd, integration of PDF from \code{-Inf} (or a point of choice) to \code{x}
#'
#' @param object An object of ecd class
#' @param x A numeric vector of \code{x}
#' @param from.x A value or a vector of starting \code{x}, default \code{-Inf}
#' @param piece.wise Logical. If \code{TRUE}, use cumulative method for large array.
#'        Default to FALSE.
#'        Use it with a scalar \code{from.x}.
#' @param f an optional extension to perform integral on function other than 1.
#'          This is for internal use only. You should use the respective wrappers.
#' @param verbose logical, display timing information, for debugging purpose.
#'
#' @return The CDF
#'
#' @keywords cdf
#'
#' @author Stephen H. Lihn
#'
#' @export
#'
#' @examples
#' d <- ecd()
#' x <- seq(-10, 10, by=1)
#' ecd.cdf(d,x)
#' ecd.cdf(d,1, from.x = -1)
### <======================================================================>
"ecd.cdf" <- function(object, x, from.x = -Inf, piece.wise = FALSE, f = NULL, verbose=FALSE)
{
    if (is.null(f)) {
    	f <- function(x) {1}
    }

    N <- length(x)
    N0 <- length(from.x)
    
    if (N0 != 1 & N0 != N) {
        stop("parameter 'from.x' must be a constant or have the same length as x!")
    }
    
    # ---------------------------------------------
    # handling piece.wise mode
    
    if (piece.wise == TRUE & N > 1) {
        if (N0 != 1) {
            stop("parameter 'from.x' must be a constant in piece.wise mode!")
        }
        if (verbose) print(paste(Sys.time(), "ecd.cdf: piece.wise"))
        
        x0 <- ecd.lag(x,1)
        x0[1] <- from.x
        pieces <- ecd.cdf(object, x, from.x=x0, f=f, verbose=verbose)
        return(cumsum(pieces))
        
    }
    
    # ---------------------------------------------
    # start the integration of x array
    
    x0 <- ecd.mpnum(object,rep(0,N)) + from.x # ensure it is an array
    
    intg <- function(z1, z2) {
        if (verbose) print(paste(Sys.time(), "ecd.cdf: intg from",
                                ecd.mp2f(z1), "to", ecd.mp2f(z2)
                         ))
        m <- integrate_pdf(object, f, z1, z2, verbose=verbose)
        if (m$message != "OK") {
            stop(paste("ecd.cdf: Failed to integrate from", ecd.mp2f(z1),
                       "to", ecd.mp2f(z2),
                       "msg:", m$message, "from ecd:", ecd.toString(object)
                ))
        }
        unname(m$value)
    }

    ecd.mpnum(object, mapply(intg, x0, x))
    
}
### <---------------------------------------------------------------------->
