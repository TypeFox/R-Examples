#' @title Black-Scholes Formula (European Option)
#' @description The famous Black-Scholes Option Pricing Formula based on the Lognormal Models. This formula can be extended to barrier options, currency options, options on futures, etc.
#' @param S The Stock Price
#' @param K The Strike Price
#' @param r The risk-free continuously compounded interest rate
#' @param T The expiration date
#' @param sigma The volatility
#' @details The Black-Scholes Formula is based on the assumption of geometric brownian motion and can be shown to satisfy the Black-Scholes Partial Differential Equation. It can be thought of as the combination of an asset-or-nothing option and a cash-or-nothing option
#' @examples blackScholes(S=41,K=40,r=0.08,T=1,sigma=0.3)
#' @export
blackScholes <- function(S, K, r, T, sigma) {  
  d1 <- (log(S/K)+(r+sigma^2/2)*T)/(sigma*sqrt(T))
  d2 <- d1 - sigma * sqrt(T)
  
  call = S*pnorm(d1) - K*exp(-r*T)*pnorm(d2)
  put  =  K*exp(-r*T) * pnorm(-d2) - S*pnorm(-d1)
  return(call)
}