#' Compute statistics of an ecd object
#' 
#' Compute statistics for m1, m2, m3, m4, mean, var, skewness, kurtosis.
#' This is used as part of ecd constructor.
#'
#' @param object an object of ecd class
#' @param asymp.q If specified, a length-one numeric as asymptotic quantile
#'     for the asymptotic statistics. There is a wrapper in \code{\link{ecd.asymp_stats}}
#' @param verbose logical, display timing information, for debugging purpose.
#'
#' @return a list of m1, m2, m3, m4, mean, var, skewness, kurtosis
#'
#' @keywords statistics
#'
#' @author Stephen H. Lihn
#'
#' @export
#'
#' @examples
#' d <- ecd(1,1)
#' ecd.stats(d)

### <======================================================================>
"ecd.stats" <- function(object, asymp.q=NULL, verbose=FALSE)
{
    m1 <- NULL
    m2 <- NULL
    m3 <- NULL
    m4 <- NULL
    
    if (is.null(asymp.q)) {
        if (! object@use.mpfr) {
            m1 <- moment(object, 1, verbose=verbose)
            m2 <- moment(object, 2, verbose=verbose)
            m3 <- moment(object, 3, verbose=verbose)
            m4 <- moment(object, 4, verbose=verbose)
        } else {
            # mpfr takes a long time, so it makes sense to use mclapply
            mf <- function(n) moment(object, n, verbose=verbose)
            mn <- c(1,2,3,4)
            mvec <- ecd.mpfr(parallel::mclapply(mn, mf))
            m1 <- mvec[1]
            m2 <- mvec[2]
            m3 <- mvec[3]
            m4 <- mvec[4]
        }
    }
    else {
        if (length(asymp.q)>1) stop("ecd.stats cannot handle vector of asymp.q")

        x1 <- qec(asymp.q, object)
        x2 <- qec(1-asymp.q, object)
        m1 <- moment(object, 1, asymp.lower=x1, asymp.upper=x2, verbose=verbose)
        m2 <- moment(object, 2, asymp.lower=x1, asymp.upper=x2, verbose=verbose)
        m3 <- moment(object, 3, asymp.lower=x1, asymp.upper=x2, verbose=verbose)
        m4 <- moment(object, 4, asymp.lower=x1, asymp.upper=x2, verbose=verbose)
    }
    
    if (abs(m1) < .Machine$double.eps) {
        m1 <- m1*0 # literally zero
    }
    if (abs(m3) < .Machine$double.eps) {
        m3 <- m3*0 # literally zero
    }
    
    var <- m2-m1^2
    
    s <- list(m1=m1, m2=m2, m3=m3, m4=m4, mean=m1, 
              var=var, stdev=sqrt(var))
             
    s$skewness <- (m3 - 3*m1*m2 + 2*m1^3 ) / var^(3/2)
    s$kurtosis <- (m4 - 4*m1*m3 + 6*m1^2*m2 - 3*m1^4 ) / var^2

    if (verbose) print(paste(Sys.time(), "ecd.stats: done, kurtosis=",
                            ecd.mp2f(s$kurtosis)
                            ))

    s
}
### <---------------------------------------------------------------------->

