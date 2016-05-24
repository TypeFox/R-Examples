#' Derives VaR of a short Black Scholes put option
#' 
#' Function derives the VaR of a Short Black Scholes put for specified 
#' confidence level and holding period, using analytical solution.
#' 
#' @param stockPrice Stock price of underlying stock
#' @param strike Strike price of the option
#' @param r Risk-free rate and is annualised
#' @param mu Mean return
#' @param sigma Volatility of the underlying stock
#' @param maturity Term to maturity and is expressed in days
#' @param cl Confidence level and is scalar
#' @param hp Holding period and is scalar and is expressed in days
#' @return Price of European put Option
#' @references Dowd, Kevin. Measuring Market Risk, Wiley, 2007.
#' 
#' Hull, John C.. Options, Futures, and Other Derivatives. 4th ed., Upper Saddle
#' River, NJ: Prentice Hall, 200, ch. 11.
#' 
#' Lyuu, Yuh-Dauh. Financial Engineering & Computation: Principles, 
#' Mathematics, Algorithms, Cambridge University Press, 2002.
#' 
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Derives VaR of a short Black Scholes put option
#'    ShortBlackScholesPutVaR(27.2, 25, .03, .12, .2, 60, .95, 40)
#'    
#' @export
ShortBlackScholesPutVaR <- function(stockPrice, strike, r, mu, sigma, 
                                   maturity, cl, hp){
  # Simplify notation
  t <- maturity/360
  hp <- hp/360
  currentOptionPrice <- BlackScholesPutPrice(stockPrice, strike, r, sigma, t)
  sStar <- exp(log(stockPrice) + (mu - (sigma ^ 2)/2) * hp - 
                 qnorm(cl, 0, 1) * sigma * sqrt(hp))
  # Critical future stock price, see Hull, p. 238
  futureOptionPriceStar <- BlackScholesPutPrice(sStar, strike, r, sigma, t - hp)
  y <- - currentOptionPrice + futureOptionPriceStar
  return(y)
}