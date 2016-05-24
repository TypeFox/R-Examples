#' Implied volatility of Black-Sholes model
#'
#' This is the standard library to calculate implied volatility \eqn{\sigma_{BS}}
#' in Black-Sholes model.
#' There is no external dependency on elliptic distribution.
#'
#' @param V numeric vector of option prices
#' @param K numeric vector of strike prices
#' @param S length-one numeric for underlying price
#' @param ttm length-one numeric for time to maturity, in the unit of days/365.
#' @param int_rate length-one numeric for risk-free rate, default to 0.
#' @param div_yield length-one numeric for dividend yield, default to 0.
#' @param otype character, specifying option type: \code{c} or \code{p}.
#' @param stop.on.na logical, to stop if fails to find solution.
#'                   Default is to use NaN and not stop.
#' @param use.mc logical, to use mclapply (default), or else just use for loop.
#'               For loop option is typically for debugging.
#'
#' @return The implied volatity \eqn{\sigma_{BS}}.
#'
#' @keywords ecop
#'
#' @export
#'
#' @examples
#' V <- c(1.8, 50)
#' K <- c(2100, 2040)
#' S <- 2089.27
#' T <- 1/365
#' y <- 0.019
#' ecop.bs_implied_volatility(V, K, S, ttm=T, div_yield=y, otype="c")
#' # expect output of 12.8886% and 29.4296%

### <======================================================================>
"ecop.bs_implied_volatility" <- function(V, K, S, ttm, int_rate=0, div_yield=0, otype="c",
                                         stop.on.na=FALSE, use.mc=TRUE)
{
    if (!(otype %in% c("c","p"))) {
        stop(paste("Unknown option type:", otype))
    }
    
    use.mpfr <- ifelse(class(V)=="mpfr" | class(K)=="mpfr", TRUE, FALSE)
    
    len = length(V)
    if (len != length(K)) {
        stop("Length of V and K must match!")
    }
    if (len > 1 & use.mc==TRUE) {
        f <- function(i) ecop.bs_implied_volatility(V[i], K[i], S, ttm,
                            int_rate, div_yield, otype=otype,
                            stop.on.na=stop.on.na)
        s1 <- parallel::mclapply(1:len, f)
        s1 <- if (use.mpfr) ecd.mpfr(s1) else simplify2array(s1)
        return(s1)
    }
    # handle length-one numeric if use.mc

    sigma = V*NaN
    for (i in 1:len)
    {
        lower = 0.000001
        upper = 4
        
        df <- function(x) {
            # x is implied volatility
            valBS <- ecop.bs_option_price(x, K[i], S, ttm, int_rate, div_yield, otype)
            valBS-V[i]
        }
        
        IV <- suppressWarnings(tryCatch(
                ecd.uniroot(df, lower = lower, upper = upper, use.mpfr = use.mpfr,
                            tol=1e-4, maxiter = 5000),
                error = function(x) NULL))
        
        if (!is.null(IV)) {
            IV$err <- ecd.mp2f(abs( df(IV$root)/V[i] ))
            if (IV$err > 1e-3 & stop.on.na) {
            	stop(paste("IV err too large, err=", IV$err, "for", otype, i,
                           "V=", ecd.mp2f(V), "K=", ecd.mp2f(K),
                           "S=", ecd.mp2f(S), "T=", ecd.mp2f(T), "\n"
                           ))
            }
        } else if (stop.on.na) {
            stop(paste("No IV found for", otype, i,
                       "V=", ecd.mp2f(V), "K=", ecd.mp2f(K),
                       "S=", ecd.mp2f(S), "T=", ecd.mp2f(T), "\n"
                       ))
        }
        
        sigma[i] <- (if(!is.null(IV)) IV$root else NaN)
    }
    sigma     
}
### <---------------------------------------------------------------------->
