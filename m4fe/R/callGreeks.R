d1 <- function(S, K, sigma, r, delta, t) {
  (log(S/K) + (r - delta + 1/2*sigma^2)*t)/(sigma*sqrt(t))
}

d2 <- function(S, K, sigma, r, delta, t) {
  d1(S, K, sigma, r, delta, t) - sigma*sqrt(t)
}

#' @title Delta (Call Greek)
#' @description A partial derivative of the Black-Scholes Equation: dc/dS (with respect to the stock price)
#' @param S The Stock Price
#' @param K The Strike Price
#' @param sigma The volatility
#' @param r The continuously compounded risk-tree interest rate
#' @param delta The annualized dividend rate
#' @param t The expiration date
#' @examples Delta_c(S=40,K=40,sigma=0.3,r=0.08,delta=0,t=91/365)
#' @export
Delta_c <- function(S, K, sigma, r, delta, t) {
  exp(-delta*t)*pnorm(d1(S, K, sigma, r, delta, t))
}

#' @title Gamma (Call Greek)
#' @description A partial derivative of the Black-Scholes Equation: dc^2/dS^2 (with respect to the stock price)
#' @param S The Stock Price
#' @param K The Strike Price
#' @param sigma The volatility
#' @param r The continuously compounded risk-tree interest rate
#' @param delta The annualized dividend rate
#' @param t The expiration date
#' @examples Gamma_c(S=40,K=40,sigma=0.3,r=0.08,delta=0,t=91/365)
#' @export
Gamma_c <- function(S, K, sigma, r, delta, t) {
  exp(-delta*t)*dnorm(d1(S, K, sigma, r, delta, t))/(S*sigma*sqrt(t))
}

#' @title Vega (Call Greek)
#' @description A partial derivative of the Black-Scholes Equation: dc/dsigma (with respect to the volatility)
#' @param S The Stock Price
#' @param K The Strike Price
#' @param sigma The volatility
#' @param r The continuously compounded risk-tree interest rate
#' @param delta The annualized dividend rate
#' @param t The expiration date
#' @examples Vega_c(S=40,K=40,sigma=0.3,r=0.08,delta=0,t=91/365)
#' @export
Vega_c <- function(S, K, sigma, r, delta, t) {
  (K*exp(-r*t)*dnorm(d1(S, K, sigma, r, delta, t))*sqrt(t))*0.01
}

#' @title Theta (Call Greek)
#' @description A partial derivative of the Black-Scholes Equation: -dc/dt (with respect to the expiration time)
#' @param S The Stock Price
#' @param K The Strike Price
#' @param sigma The volatility
#' @param r The continuously compounded risk-tree interest rate
#' @param delta The annualized dividend rate
#' @param t The expiration date
#' @examples Theta_c(S=40,K=40,sigma=0.3,r=0.08,delta=0,t=91/365)
#' @export
Theta_c <- function(S, K, sigma, r, delta, t) {
  (delta*S*exp(-delta*t)*pnorm(d1(S, K, sigma, r, delta, t)) - K*exp(-r*t)*r*pnorm(d2(S, K, sigma, r, delta, t)) - K*exp(-r*t)*sigma/(2*sqrt(t))*dnorm(d2(S, K, sigma, r, delta, t)))/365
}

#' @title Rho (Call Greek)
#' @description A partial derivative of the Black-Scholes Equation: dc/dr (with respect to the interest rate)
#' @param S The Stock Price
#' @param K The Strike Price
#' @param sigma The volatility
#' @param r The continuously compounded risk-tree interest rate
#' @param delta The annualized dividend rate
#' @param t The expiration date
#' @examples Rho_c(S=40,K=40,sigma=0.3,r=0.08,delta=0,t=91/365)
#' @export
Rho_c <- function(S, K, sigma, r, delta, t) {
  (t*K*exp(-r*t)*pnorm(d2(S, K, sigma, r, delta, t)))*0.01
}

#' @title Psi (Call Greek)
#' @description A partial derivative of the Black-SCholes Equation: dc/ddelta (with respect to the dividend rate)
#' @param S The Stock Price
#' @param K The Strike Price
#' @param sigma The volatility
#' @param r The continuously compounded risk-tree interest rate
#' @param delta The annualized dividend rate
#' @param t The expiration date
#' @examples Psi_c(S=40,K=40,sigma=0.3,r=0.08,delta=0,t=91/365)
#' @export
Psi_c <- function(S, K, sigma, r, delta, t) {
  (-t*S*exp(-delta*t)*pnorm(d1(S, K, sigma, r, delta, t)))*0.01
}