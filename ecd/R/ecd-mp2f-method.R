#' Wrapper to convert mpfr to numeric
#'
#' Convert mpfr to numeric primarily for display messages.
#'
#' @param x an object of \code{mpfr} class. If \code{x} is numeric class,
#'          it will be passed through.
#'
#' @return a numeric vector
#'
#' @keywords utility
#'
#' @export
#'
#' @importFrom Rmpfr asNumeric
#'
#' @examples
#' x <- ecd.mp2f(ecd.mpfr(c(1,2,3)))
### <======================================================================>
"ecd.mp2f" <- function(x)
{
    c <- class(x)
    if (c=="numeric") return(x)
    if (c=="integer") return(as.numeric(x))
    if (c=="mpfr") {
        n <- length(x)
        z <- rep(NaN,n)
        for (i in 1:n) {
            z[i] <- Rmpfr::asNumeric(x[i])
        }
        return(z)
    }
    stop(paste("Unsupported type for mp2f:", c))
}
### <---------------------------------------------------------------------->
