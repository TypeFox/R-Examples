#' Utility to shift a vector of numeric or mpfr
#'
#' This utility is basically the same as \code{Hmisc::Lag},
#' but it handles mpfr vector properly.
#'
#' @param x a vector of numeric or mpfr
#' @param shift integer, cells to shift
#' @param na.omit logical, whether to remove the NAs
#'
#' @return the shifted vector
#'
#' @keywords utility
#'
#' @export
#'
#' @examples
#' x <- ecd.lag(c(1,2,3))
#' y <- ecd.lag(ecd.mpfr(c(1,2,3)))
### <======================================================================>
ecd.lag <- function(x, shift=1, na.omit=FALSE)
{
    ## Lags vector x shift observations, padding with NAs
    ## preserving attributes of x
    
    xLen <- length(x)
    if (shift == 0) return(x)
    
    ret <- rep(NaN, xLen)
    if (class(x) == "mpfr") {
        ret <- ecd.mpfr(ret)
    }    
    if (abs(shift) < xLen) {
        if (shift > 0) ret[-(1:shift)] <- x[1:(xLen - shift)]
        else ret[1:(xLen+shift)] <- x[(1-shift):xLen]
    }
    if (na.omit) {
        ret <- ret[!is.na(ret)]
    }
    return(ret)
}
### <---------------------------------------------------------------------->
