#' @title European Call Option Price (Recursive)
#' @description Determines the price of a European Call Option by using a recursive relationship in the binomial tree
#' @param S The Stock Price
#' @param K The Strike Price
#' @param sigma The volatility
#' @param r The risk-free continuously compounded interest rate
#' @param delta The annualized dividend rate
#' @param h the number of periods betwen 0 and T, where each period is of length 1/h
#' @param T the expiration time
#' @details Uses formulas for u and d presented by Cox in his forward tree model. Note that the option price converges to the Black-Scholes option price
#' @examples recursiveTree(S=41,K=40,sigma=0.3,r=0.08,delta=0,h=5,T=1)
#' @export
recursiveTree <- function(S, K, sigma, r, delta, h, T=1) {
  if(T!=1) {  r = r*T; h = h*T; delta = delta*T; sigma = sigma*sqrt(T) }
  
  u = exp((r-delta)/h + sigma*sqrt(1/h)); d = exp((r-delta)/h - sigma*sqrt(1/h))
  p = (exp((r-delta)/h) - d)/(u-d)
  
  cost <- function(node) {
    if(length(node)==h) return(max(0,Sundn(node) - K))
    else return(max(exp(-(r-delta)/h)*(p*cost(c(node,1)) + (1-p)*cost(c(node,0))),0))
  }
  
  Sundn <- function(node) {
    un = sum(node); dn = length(node) - un
    S*u^un*d^dn
  }
  
  cost(numeric(0))
}