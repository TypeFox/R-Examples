#'  Geometric average asian options
#'
#' @description Pricing functions for European Asian options based on
#'     geometric averages. \code{geomavgprice} and
#'     \code{geomavgstrike} compute analytical prices of geometric
#'     Asian options using the modified Black-Scholes
#'     formula. 
#'

#' @title Geometric average price
#' @export
#' @family Asian
#' @name geomavgprice
#' @param s Stock price
#' @param k Strike price of the option. In the case of average strike
#'     options, \code{k/s} is the multiplier for the average
#' @param v Volatility of the stock, defined as the annualized
#'     standard deviation of the continuously-compounded return
#' @param r Annual continuously-compounded risk-free interest rate
#' @param tt Time to maturity in years
#' @param d Dividend yield, annualized, continuously-compounded
#' @param m Number of prices in the average calculation
#' @param cont Boolean which when TRUE denotes continuous averaging
#' 
#' @examples
#' s=40; k=40; v=0.30; r=0.08; tt=0.25; d=0; m=3;
#' geomavgprice(s, k, v, r, tt, d, m)
#' 
#' @description Prices of geometric average-price call and put options
#' @usage
#' geomavgprice(s, k, v, r, tt, d, m, cont=FALSE)
#' @return Call and put prices as a vector
#'
geomavgprice <- function(s, k, v, r, tt, d, m, cont=FALSE) {
    if (cont) {
        siga <- v / (3^0.5)
        da <- 0.5 * (r + d + v^2 / 6)
    } else {
        siga <- v * ((m + 1) * ((2 * m) + 1) / 6)^0.5 / m 
        da <- 0.5 * (r * (m - 1) / m + (d + 0.5 * v^2) *
                     (m + 1) / m - (v / m)^2 *
                    (m + 1) * (2 * m + 1) / 6)
    }
    return(c(Call=bscall(s, k, siga, r, tt, da),
             Put=bsput(s, k, siga, r, tt, da)))
}

#' @export
#' @family Asian
#' @name geomavgstrike
#' @description Prices of geometric average-strike call and put options
#' @param km The strike mutiplier, relative to the initial stock
#'     price, for an average price payoff. If the initial stock price
#'     is \code{s = 120} and \code{km = 115}, the payoff for an
#'     average strike call is \deqn{Payoff = max(ST - km/s*SAvg, 0)}.
#' @usage
#' geomavgstrike(s, km, v, r, tt, d, m, cont=FALSE)
#' @examples
#' s=40; km=40; v=0.30; r=0.08; tt=0.25; d=0; m=3;
#' @return Vector of call and put prices for geometric average strike
#'     options
#' geomavgstrike(s, km, v, r, tt, d, m)
#' @title Geometric average-strike options
#' @inheritParams geomavgprice
geomavgstrike <- function(s, km, v, r, tt, d, m, cont=FALSE) {
    if (cont) {
        siga <- v / (3 ^ 0.5)
        da <- 0.5 * (r + d + v ^ 2 / 6)
        rho <- 0.5 * (3 ^ 0.5)
    } else {
        siga <- v * ((m + 1) * (2 * m + 1) / 6) ^ 0.5 / m
        da <- 0.5 * (r * (m - 1) / m + (d + 0.5 * v ^ 2) * (m + 1) / m -
                      (v / m) ^ 2 * (m + 1) * (2 * m + 1) / 6)
        rho <- 0.5 * (6 * (m + 1) / (2 * m + 1)) ^ 0.5
    }
    vol <-  (siga ^ 2 + v ^ 2 - 2 * rho * siga * v) ^ 0.5
    return(c(Call=bscall(s, km, vol, da, tt, d),
           Put=bsput(s, km, vol, da, tt, d)))
}

