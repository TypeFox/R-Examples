#' @title European Option Pricing Details
#' @description Returns a tree structure with the D (Delta: amount of shares to buy), B (lend: borrowing amount), S (stock price), C (option value) at each node
#' @param n the number of periods
#' @param call TRUE if call option, FALSE if put option
#' @details Uses a list structure to store the data. If n is large, the object is very large and takes up a lot of memory
#' @examples treeDetails(n=10, call=TRUE)
#' @export
treeDetails <- function(n=2, call=TRUE) {
  S = 41; K = 40; sigma = 0.3; r = 0.08; delta = 0; h = n; T = 1
  if(T!=1) r = r*T; h = h*T; delta = delta*T; sigma = sigma*sqrt(T)
  
  type=call
  factor = ifelse(type==TRUE, -1, 1)
  
  u = exp((r-delta)/h + sigma*sqrt(1/h))
  d = exp((r-delta)/h - sigma*sqrt(1/h))  
  p = (exp((r-delta)/h) - d)/(u-d)
  
  tree = list()
  for(i in h:0) {
    node <- numeric(0)
    value <- numeric(0)
    Delta <- numeric(0)
    lend <- numeric(0)
    
    for(j in 0:i) {
      node = c(node, S*u^(i-j)*d^j)
      if(i==h) {
        value = c(value, max(factor*(K-S*u^(i-j)*d^j), 0))
        Delta = NA
        lend = NA
      } else {
        Delta = c(Delta, (tree[[i+2]]$C[j+1] - tree[[i+2]]$C[j+2])/(tree[[i+2]]$S[j+1] - tree[[i+2]]$S[j+2]))
        lend = c(lend, exp(-(r-delta)/h)*(u*tree[[i+2]]$C[j+2] - d*tree[[i+2]]$C[j+1])/(u-d))
        value = c(value, max(exp(-(r-delta)/h)*(p*tree[[i+2]]$C[j+1] + (1-p)*tree[[i+2]]$C[j+2]),0))
      }
    }
    level = list(S=node, C=value, D=Delta, B=lend)
    tree[[i+1]] = level
  }
  return(tree)
}