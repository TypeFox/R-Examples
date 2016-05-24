#' Price of European Call Option
#' 
#' Derives the price of European call option using the Black-Scholes approach
#' 
#' @param stockPrice Stock price of underlying stock
#' @param strike Strike price of the option
#' @param rf Risk-free rate and is annualised
#' @param sigma Volatility of the underlying stock
#' @param t The term to maturity of the option in years
#' @return Price of European Call Option
#' @references Dowd, Kevin. Measuring Market Risk, Wiley, 2007.
#' 
#' Hull, John C.. Options, Futures, and Other Derivatives. 5th ed., p. 246.
#' 
#' Lyuu, Yuh-Dauh. Financial Engineering & Computation: Principles, 
#' Mathematics, Algorithms, Cambridge University Press, 2002.
#' 
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Estimates the price of an American Put
#'    BlackScholesCallPrice(27.2, 25, .03, .2, 60)
#'    
#' @export
BlackScholesCallPrice <- function(stockPrice, strike, rf, sigma, t){
  S <- stockPrice
  X <- strike
  # d terms
  d1 <- (log(S/X) + (rf + (sigma ^ 2) / 2) * t) / (sigma * sqrt(t))
  d2 <- d1 - sigma * sqrt(t)
  # Option price
  y <- stockPrice * pnorm(d1) - exp(- rf * t) * strike * pnorm(d2)
  return(y)
}