#' Option generating function of ecd
#'
#' Option generating function (OGF) of ecd.
#' For call, it is integration of \eqn{(e^z-e^k) P(z)} for \code{z} from \code{k} to \code{Inf}.
#' For put, it is integration of \eqn{(e^k-e^z) P(z)} for \code{z} from \code{-Inf} to \code{k}.
#'
#' @param object an object of ecd class
#' @param k a numeric vector of log-strike
#' @param otype character, specifying option type: \code{c} or \code{p}.
#' @param unit.sigma logical, transforming to unit sigma to achieve greater stability.
#' @param verbose logical, display timing information, for debugging purpose.
#'
#' @return The option price normalized by underlying
#'
#' @keywords option
#'
#' @author Stephen H. Lihn
#'
#' @export
#'
#' @examples
#' d <- ecd(0, 0, sigma=0.01)
#' k <- seq(-0.1, 0.1, by=0.01)
#' ecd.ogf(d, k, "c")
### <======================================================================>
"ecd.ogf" <- function(object, k, otype="c", unit.sigma=FALSE, verbose=FALSE)
{
    # TODO this can become unstable when the result is the substraction of two large numbers
    stopifnot(is.numericMpfr(k))

    if (verbose) print(paste(Sys.time(), "ecd.ogf: start",
                            "k=", ecd.mp2f(k),
                            "otype=", otype
                      ))

    if (otype=="p") {
        fn_exp <- function(z) exp(z)
        imgf2 <- ecd.cdf(object, k, f = fn_exp, piece.wise=TRUE, verbose=verbose)
        cdf <- ecd.cdf(object, k, piece.wise=TRUE, verbose=verbose)
        
        if (verbose) print(paste(Sys.time(), "ecd.ogf: done"))
        
        return(exp(k)*cdf - imgf2)
    }
    
    if (otype=="c") {
        imgf <- ecd.imgf(object, k, unit.sigma= unit.sigma, verbose=verbose)
        ccdf <- ecd.ccdf(object, k, piece.wise=TRUE, verbose=verbose)

        if (verbose) print(paste(Sys.time(), "ecd.ogf: done"))
        
        return(imgf - exp(k)*ccdf)
    }
    
    stop(paste("Unknown option type:", otype))

}
### <---------------------------------------------------------------------->
