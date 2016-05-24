#' @title European Option Price (Basic Binomial Tree Formula)
#' @description Uses Cox's Formula for finding the price of a European call option.
#' @param n the number of periods
#' @param call TRUE if call option, FALSE if put option
#' @details This formula assumes that nodes do overlap and path dependency is not needed as in an Asian option.
#' @examples treeBasic(n=10, call=TRUE)
#' @export
treeBasic <- function(n=3, call=TRUE) {
  S = 41; K = 40; sigma = 0.3; r = 0.08; delta = 0; h = n; T = 1
  if(T!=1) r = r*T; h = h*T; delta = delta*T; sigma = sigma*sqrt(T)
  
  type=call
  factor = ifelse(type==TRUE, -1, 1)
  
  u = exp((r-delta)/h + sigma*sqrt(1/h))
  d = exp((r-delta)/h - sigma*sqrt(1/h))  
  p = (exp((r-delta)/h) - d)/(u-d)
  
  for(i in h:0) {
    value <- numeric(0)
    for(j in 0:i) {
      if(i==h) {
        value = c(value, max(factor*(K-S*u^(i-j)*d^j), 0))
      } else {
        value = c(value, max(exp(-(r-delta)/h)*(p*last[j+1] + (1-p)*last[j+2]),0))
      }
    }
    last = value
  }
  return(last)
}