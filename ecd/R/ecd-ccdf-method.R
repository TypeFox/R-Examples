#' Complementary CDF of ecd
#'
#' Complementary CDF of ecd, integration of PDF from \code{x} to \code{Inf}
#'
#' @param object An object of ecd class
#' @param x A numeric vector of \code{x}
#' @param to.x A value or a vector of starting \code{x}, default \code{Inf}
#'             This is for internal use only.
#' @param piece.wise Logical. If \code{TRUE}, use cumulative method for large array.
#'        Default to FALSE.
#'        Use it with a scalar \code{to.x}.
#' @param f an optional extension to perform integral on function other than 1.
#'          This is for internal use only. You should use the respective wrappers.
#' @param verbose logical, display timing information, for debugging purpose.
#'
#' @return The CCDF
#'
#' @keywords cdf
#'
#' @author Stephen H. Lihn
#'
#' @export
#'
#' @examples
#' d <- ecd()
#' x <- seq(0, 10, by=1)
#' ecd.ccdf(d,x)
### <======================================================================>
"ecd.ccdf" <- function(object, x, to.x = Inf, piece.wise = FALSE, f = NULL, verbose=FALSE)
{
    if (is.null(f)) {
    	f <- function(x) {1}
    }
    
    N <- length(x)
    N0 <- length(to.x)
    
    if (N0 != 1 & N0 != N) {
        stop("parameter 'to.x' must be a constant or have the same length as x!")
    }
    
    # ---------------------------------------------
    # handling piece.wise mode
    
    if (piece.wise == TRUE & N > 1) {
        if (N0 != 1) {
            stop("parameter 'to.x' must be a constant in piece.wise mode!")
        }
        if (verbose) print(paste(Sys.time(), "ecd.ccdf: piece.wise"))
       
        xr <- rev(x)
        x0 <- ecd.lag(xr,1)
        x0[1] <- to.x
        pieces <- ecd.ccdf(object, xr, to.x=x0, f=f, verbose=verbose)
        return(rev(cumsum(pieces)))
        
    }
    
    # ---------------------------------------------
    # start the integration of x array
    
    x0 <- ecd.mpnum(object,rep(0,N)) + to.x # ensure it is an array

    intg2 <- function(z1, z2) {
        
        if (verbose) print(paste(Sys.time(), "ecd.ccdf: intg2 from",
                                ecd.mp2f(z1), "to", ecd.mp2f(z2)
                         ))

        m <- integrate_pdf(object, f, z1, z2, verbose=verbose)
        if (m$message != "OK") {
            stop(paste("ecd.ccdf: Failed to integrate from", ecd.mp2f(z1),
                       "to", ecd.mp2f(z2),
                       "msg:", m$message, "from ecd:", ecd.toString(object)
                ))
        }
        unname(m$value)
    }
    
    ecd.mpnum(object, mapply(intg2, x, x0))

}
### <---------------------------------------------------------------------->
