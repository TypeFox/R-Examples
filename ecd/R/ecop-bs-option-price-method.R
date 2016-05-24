#' Calculate option price from implied volatility in Black-Sholes model
#'
#' This is the standard library to calculate option price from
#' implied volatility \eqn{\sigma_{BS}} in Black-Sholes model.
#' There is no external dependency on elliptic distribution.
#'
#' @param ivol numeric vector of implied volatility
#' @param K numeric vector of strike prices
#' @param S length-one numeric for underlying price
#' @param ttm length-one numeric for time to maturity, in the unit of days/365.
#' @param int_rate length-one numeric for risk-free rate, default to 0.
#' @param div_yield length-one numeric for dividend yield, default to 0.
#' @param otype character, \code{c} or \code{p}. Default is \code{c}.
#'
#' @return The call/put prices
#'
#' @keywords ecop
#'
#' @export ecop.bs_option_price
#' @export ecop.bs_call_price
#' @export ecop.bs_put_price
#'
#' @importFrom Rmpfr pnorm
#'
#' @examples
#' ivol <- c(0.128886, 0.294296) 
#' K <- c(2100, 2040)
#' S <- 2089.27
#' T <- 1/365
#' y <- 0.019
#' ecop.bs_option_price(ivol, K, S, ttm=T, div_yield=y, otype="c")
#' # expect output of c(1.8, 50)

### <======================================================================>
ecop.bs_option_price <- function(ivol, K, S, ttm, int_rate=0, div_yield=0, otype="c")
{
    if (!(otype %in% c("c","p"))) {
        stop(paste("Unknown option type:", otype))
    }

    T = ttm
    r = int_rate
    
    d1 <- (log(S/K) + (r - div_yield + (ivol^2)/2) * T) / (ivol * sqrt(T))
    d2 <- d1 - ivol * sqrt(T)
    use.mpfr <- ifelse(class(d1)=="mpfr" | class(d2)=="mpfr", TRUE, FALSE)

    sgn <- if (otype=="p") -1 else 1
    d1 <- ecd.mp1* sgn* d1
    d2 <- ecd.mp1* sgn* d2

    SD <- S * exp(-div_yield*T)
    KD <- K * exp(-r*T)

    # Rmpfr pnorm doesn't like NaN
    # use MPFR pnorm to enhance accuracy
    pnorm2 <- function(d) {
        p <- rep(NaN, length(d))*ecd.mp1
        di <- which(!is.na(d))
        p[di] <- Rmpfr::pnorm(d[di])
        p
    }

    V <- sgn* (SD * pnorm2(d1) - KD * pnorm2(d2))
    if (use.mpfr) return(V) else return(ecd.mp2f(V))
    
}
### <---------------------------------------------------------------------->
#' @rdname ecop.bs_option_price
ecop.bs_call_price <- function(ivol, K, S, ttm, int_rate=0, div_yield=0)
{
    ecop.bs_option_price(ivol, K, S, ttm, int_rate, div_yield, otype="c")
}
#' @rdname ecop.bs_option_price
ecop.bs_put_price <- function(ivol, K, S, ttm, int_rate=0, div_yield=0)
{
    ecop.bs_option_price(ivol, K, S, ttm, int_rate, div_yield, otype="p")
}






