#' Compute the moment of ecd via integration
#' 
#' Compute the moment of ecd via integration between \code{-Inf} and \code{Inf}.
#' The \code{asymp.lower} and \code{asymp.upper} parameters are used for
#' asymptotic statistics, to study the effect of finite observations.
#'
#' @method moment ecd
#'
#' @param object an object of ecd class
#' @param order numeric. Order of the moment to be computed
#' @param center logical. If set to \code{TRUE}, calculate central moments.
#'               Default: \code{FALSE}.
#' @param asymp.lower numeric, lower bound for asymptotic statistics, default: \code{-Inf}.
#' @param asymp.upper numeric, upper bound for asymptotic statistics, default: \code{Inf}.
#' @param verbose logical, display timing information, for debugging purpose.
#'
#' @return Numeric. The moment.
#'
#' @keywords moment
#'
#' @author Stephen H. Lihn
#'
#' @export moment
#'
#' @examples
#' d <- ecd()
#' moment(d, 2)


### <======================================================================>
"moment.ecd" <- function(object, order, center=FALSE,
                         asymp.lower=-Inf, asymp.upper=Inf, verbose=FALSE)
{
    if (verbose) print(paste(Sys.time(), "moment.ecd: order", order, "center", center))

    # analytic result: odd moment = mu if beta = 0
    #if (object@beta == 0 & floor(order/2)*2 != order) {
    #    return(ecd.mpnum(object, ifelse(center, 0, object@mu)))
    #}
    
    # handle central moments
    c <- ifelse(center, moment(object, 1, verbose=verbose), 0)
    mnt <- function(x){(x-c)^order}
    
    # even moment, do it directly
    if (floor(order/2)*2==order) {
        m <- integrate_pdf(object, mnt, asymp.lower, asymp.upper, verbose=verbose)
        if (m$message != "OK") {
            stop(paste("moment.ecd: Failed to integrate order", order,
                       "msg:", m$message, "from ecd:", ecd.toString(object)
            ))
        }
        unname(m$value)
    } else {
        # odd moment, split between mu
        mu <- object@mu
        m1 <- integrate_pdf(object, mnt, asymp.lower, mu, verbose=verbose)
        m2 <- integrate_pdf(object, mnt, mu, asymp.upper, verbose=verbose)

        if (m1$message != "OK") {
            stop(paste("moment.ecd: Failed to integrate lower-mu for order", order,
                       "msg:", m1$message, "from ecd:", ecd.toString(object)
            ))
        }

        if (m2$message != "OK") {
            stop(paste("moment.ecd: Failed to integrate mu-upper for order", order,
                       "msg:", m2$message, "from ecd:", ecd.toString(object)
            ))
        }

        unname(m1$value+m2$value)
    }
}
### <---------------------------------------------------------------------->
#' @rdname moment.ecd
setGeneric("moment", function(object, order, center=FALSE,
    asymp.lower=-Inf, asymp.upper=Inf, verbose=FALSE) standardGeneric("moment"))
#' @rdname moment.ecd
setMethod("moment", signature("ecd"), moment.ecd)
### <---------------------------------------------------------------------->
