#' @title Black-Scholes implied volatility and price
#'
#' @description \code{bscallimpvol} and \code{bsputimpvol} compute
#' Black-Scholes implied volatilties. The functions \code{bscallimps}
#' and \code{bsputimps}, compute stock prices implied by a given
#' option price, volatility and option characteristics.
#'
#' @name implied
#' @aliases bscallimpvol bsputimpvol bscallimps bsputimps
#'
#' @return Implied volatility (for the "impvol" functions) or implied
#' stock price (for the "impS") functions.
#'
#' @importFrom stats uniroot
#' 
#' @usage
#' bscallimpvol(s, k, r, tt, d, price)
#' bsputimpvol(s, k, r, tt, d, price)
#' bscallimps(s, k, v, r, tt, d, price)
#' bsputimps(s, k, v, r, tt, d, price)
#' 
#'
#' @param s Stock price
#' @param k Strike price of the option
#' @param v Volatility of the stock, defined as the annualized
#' standard deviation of the continuously-compounded return
#' @param r Annual continuously-compounded risk-free interest rate
#' @param tt Time to maturity in years
#' @param d Dividend yield, annualized, continuously-compounded
#' @param price Option price when computing an implied value
#' 
#' @details Returns a scalar or vector of option prices, depending on
#' the inputs
#'
#' @note Implied volatilties and stock prices do not exist if the
#' price of the option exceeds no-arbitrage bounds. For example, if
#' the interest rate is non-negative, a 40 strike put cannot have a
#' price exceeding $40.
#'
#' @examples
#' s=40; k=40; v=0.30; r=0.08; tt=0.25; d=0;
#' bscallimpvol(s, k, r, tt, d, 4)
#' bsputimpvol(s, k, r, tt, d, 4)
#' bscallimps(s, k, v, r, tt, d, 4)
#' bsputimps(s, k, v, r, tt, d, 4)
#' 

.tol <- .Machine$double.eps^0.5

#' @export
bscallimpvol <- function(s, k, r, tt, d, price) {
    ## this function is not vectorized
    if (price <= s*exp(-d*tt)-k*exp(-r*tt)) {
        print("Option price violates minimum bound")
        break
    } else if (price > s*exp(-d*tt)) {
        print("Option price violates maximum bound")
    } else {
        f <- function(v, s, k, r, tt, d, price) {
            return(bscall(s, k, v, r, tt, d) - price)
        }
        x <- uniroot(f, c(0.001,1000), s, k, r, tt, d, price, tol=.tol)
        return(x$root)
    }
}

## add tolerance!
#' @export
bsputimpvol <- function(s, k, r, tt, d, price) {
    ## this function is not vectorized
    x <- bscallimpvol(s, k, r, tt, d,
                      price + s*exp(-d*tt) - k*exp(-r*tt))
    return(x)
}

#' @export
bscallimps <- function(s, k, v, r, tt, d, price) {
    f <- function(s, k, v, r, tt, d, price) {
        return(bscall(s, k, v, r, tt, d) - price)
    }
    x <- uniroot(f, c(0.001,10000), k, v, r, tt, d, price, tol=.tol)
    return(x$root)
}

#' @export
bsputimps <- function(s, k, v, r, tt, d, price) {
    f <- function(s, k, v, r, tt, d, price) {
        return(bsput(s, k, v, r, tt, d) - price)
    }
    x <- uniroot(f, c(0.001,10000), k, v, r, tt, d, price, tol=.tol)
    return(x$root)
}

