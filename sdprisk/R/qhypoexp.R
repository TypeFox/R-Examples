#' p-quantile from the hypo-exponential distribution
#' 
#' The p-quantile from the hypo-exponential distribution is computed implicitly by
#' computing the root of H(x) = F(x) - p, where F is the distribution function of the
#' hypo-exponential distribution
#' 
#' @param p vector of probabilities
#' @param rate vector of (unique) rates 
#' @param interval search interval in which \code{uniroot} searches for the quantile
#' @return p-quantile from the hypo-exponential distribution
#' @author Sebastian Szugat \email{szugat@@statistik.tu-dortmund.de}

qhypoexp <- function(p, rate, interval = c(0.0, 1.0e+10)) {
    stopifnot(all(is.finite(rate)))

    invalid.p <- (p < 0 | p > 1)

    if (any(invalid.p)) {
        warning(paste('Some elements of', sQuote('p'), 'are outside of the unit interval.'))
        is.na(p) <- invalid.p
    }

    H <- function(x, p) {
        phypoexp(q = x, rate = rate) - p
    }

    auxfun <- function(p) {
        tryCatch(expr  = uniroot(f         = H,
                                 interval  = interval,
                                 p         = p,
                                 extendInt = 'upX',
                                 tol       = .Machine$double.eps^0.5)$root,
                 error = function(err) NA_real_)
    }

    vapply(X         = p,
           FUN       = auxfun,
           FUN.VALUE = numeric(1L))

}
