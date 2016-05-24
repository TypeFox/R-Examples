#' Compute asymptotic statistics of an ecd object
#'
#' The main API for asymptotic statistics.
#' It follows the same definition of moments, except the integral of PDF
#' is limited to a range of quantile. That is to truncate the tails.
#' The asymptotic kurtosis is also called truncated kurtosis.
#'
#' @param object an object of ecd class with quantile
#' @param q numeric vector of quantiles
#'
#' @return a list of stats list, or a vector of kurtosis
#'
#' @keywords statistics
#'
#' @export ecd.asymp_stats
#' @export ecd.asymp_kurtosis
#'
#' @examples
#' \dontrun{
#'     d <- ecd(1,1, with.quantile=TRUE)
#'     q <- 0.01
#'     ecd.asymp_stats(d,q)
#'     ecd.asymp_kurtosis(d,q)
#' }
### <======================================================================>
"ecd.asymp_stats" <- function(object, q)
{
    f <- function(q) ecd.stats(object, asymp.q=q)
    lapply(q,f)
}
### <---------------------------------------------------------------------->
#' @rdname ecd.asymp_stats
"ecd.asymp_kurtosis" <- function(object, q) {
    s <- ecd.asymp_stats(object,q)
    f <- function(s) s$kurtosis
    sapply(s,f)
}

