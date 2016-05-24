#' @title Black-Scholes option pricing
#'
#' @description \code{bscall} and \code{bsput} compute Black-Scholes
#' call and put prices. The functions \code{assetcall},
#' \code{assetput}, \code{cashcall}, and \code{cashput} provide the
#' prices of binary options that pay a share (the asset options) or $1
#' (the cash options) if at expiration the asset price exceeds the
#' strike (the calls) or is below the strike (the puts). We have the
#' identities
#'
#' \code{bscall(s, k, v, r, tt, d)
#'   = assetcall(s, k, v, r, tt, d) - k*cashcall(s, k, v, r, tt, d)}
#'
#'
#' @name blksch
#' @aliases bscall bsput assetcall assetput cashcall cashput 
#' @importFrom stats pnorm
#' @return A Black-Scholes option price. If more than one argument is a
#' vector, the recycling rule determines the handling of the inputs
#'
#' @usage
#' bscall(s, k, v, r, tt, d)
#' bsput(s, k, v, r, tt, d)
#' assetcall(s, k, v, r, tt, d)
#' cashcall(s, k, v, r, tt, d)
#' assetput(s, k, v, r, tt, d)
#' cashput(s, k, v, r, tt, d)
#' 
#'
#' @param s Stock price
#' @param k Strike price of the option
#' @param v Volatility of the stock, defined as the annualized
#' standard deviation of the continuously-compounded return
#' @param r Annual continuously-compounded risk-free interest rate
#' @param tt Time to maturity in years
#' @param d Dividend yield, annualized, continuously-compounded
#' 
#' @details Returns a scalar or vector of option prices, depending on
#' the inputs
#'
#' @note It is possible to specify the inputs either in terms of an
#' interest rate and a "dividend yield" or an interest rate and a
#' "cost of carry". In this package, the dividend yield should be
#' thought of as the cash dividend received by the owner of the
#' underlying asset, \emph{or} (equivalently) as the payment received
#' if the owner were to lend the asset.
#'
#' There are other option pricing packages available for R, and these
#' may use different conventions for specifying inputs. In fOptions,
#' the dividend yield is replaced by the generalized cost of carry,
#' which is the net payment required to fund a position in the
#' underlying asset. If the interest rate is 10\% and the dividend
#' yield is 3\%, the generalized cost of carry is 7\% (the part of the
#' interest payment not funded by the dividend payment). Thus, using
#' the \code{GBS} function from fOptions, these two expressions return
#' the same price:
#'
#' \code{bscall(s, k, v, r, tt, d)}
#'
#' \code{fOptions::GBSOption('c', S=s, K=k, Time=tt, r=r, b=r-d, sigma=v) }
#'
#'
#' @examples
#' s=40; k=40; v=0.30; r=0.08; tt=0.25; d=0;
#' bscall(s, k, v, r, tt, d)
#'
#' ## following returns the same price as previous
#' assetcall(s, k, v, r, tt, d) - k*cashcall(s, k, v, r, tt, d)
#'
#' ## return option prices for different strikes
#' bsput(s, k=38:42, v, r, tt, d)

#' @export
bscall <- function(s, k, v, r, tt, d) {
    callp <- s*exp(-d*tt)*.nd1(s, k, v, r, tt, d) -
        k*exp(-r*tt)*.nd2(s, k, v, r, tt, d)
    return(callp)
}

#' @export
bsput <- function(s, k, v, r, tt, d) {
    putp <- bscall(s, k, v, r, tt, d) + k*exp(-r*tt) -
        s*exp(-d*tt)
    return(putp)
}

#' @export
assetcall <- function(s, k, v, r, tt, d) {
    price <- s*exp(-d*tt)*pnorm(.d1(s, k, v, r, tt, d))
    return(price)
}

#' @export
cashcall <- function(s, k, v, r, tt, d) {
    price <- exp(-r*tt)*pnorm(.d2(s, k, v, r, tt, d))
    return(price)
}

#' @export
assetput <- function(s, k, v, r, tt, d) {
    price <- s*exp(-d*tt)*pnorm(-.d1(s, k, v, r, tt, d))
    return(price)
}

#' @export
cashput <- function(s, k, v, r, tt, d) {
    price <- exp(-r*tt)*pnorm(-.d2(s, k, v, r, tt, d))
    return(price)
}


.d1 <- function(s, k, v, r, tt, d)
    (log(s/k) + (r-d+v^2/2)*tt)/(v*sqrt(tt))

.d2 <- function(s, k, v, r, tt, d)
    .d1(s, k, v, r, tt, d) - v*tt^(0.5)

.nd1 <- function(s, k, v, r, tt, d)
    pnorm(.d1(s, k, v, r, tt, d))

.nd2 <- function(s, k, v, r, tt, d)
    pnorm(.d2(s, k, v, r, tt, d))

