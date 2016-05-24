#' @title Quotient option valuation via Black-Scholes (BS) model
#' @description Quotient Option via Black-Scholes (BS) model
#' @author Chengwei Ge, Department of Statistics, Rice University, Spring 2015
#' 
#' @param o An object of class \code{OptPx}
#' @param I1 A spot price of the underlying security 1 (usually I1)
#' @param I2 A spot price of the underlying security 2 (usually I2)
#' @param g1 Payout rate of the first stock
#' @param g2 Payout rate of the 2nd stock
#' @param sigma1 a vector of implied volatilities for the associated security 1
#' @param sigma2 a vector of implied volatilities for the associated security 2 
#' @param rho is the correlation between asset 1 and asset 2
#' @return A list of class \code{QuotientBS} consisting of the original \code{OptPx} object 
#' and the option pricing parameters \code{I1},\code{I2}, \code{Type}, \code{isForeign}, and \code{isDomestic}
#' as well as the computed price \code{PxBS}.
#' 
#' @references Zhang Peter G., \emph{Exotic Options}, 2nd, 1998. \url{http://amzn.com/9810235216}.
#' @examples
#' (o = QuotientBS())$PxBS
#' 
#' o = OptPx(Opt(Style = 'Quotient', Right = "Put"), r= 0.05)
#' (o = QuotientBS(o, I1=100, I2=100, g1=0.04, g2=0.03, sigma1=0.18,sigma2=0.15, rho=0.75))$PxBS
#' 
#' o = OptPx(Opt(Style = 'Quotient',  Right = "Put", ttm=1, K=1), r= 0.05)
#' QuotientBS(o, I1=100, I2=100, g1=0.04, g2=0.03, sigma1=0.18,sigma2=0.15, rho=0.75) 
#' 
#' o = OptPx(Opt(Style = 'Quotient',  Right = "Call", ttm=1, K=1), r= 0.05)
#' QuotientBS(o, I1=100, I2=100, g1=0.04, g2=0.03, sigma1=0.18,sigma2=0.15, rho=0.75) 
#' @export
#' 
QuotientBS=function(o=OptPx(Opt(Style='Quotient')),I1=100, I2=100, g1=0.04,g2=0.03,sigma1=0.18,sigma2=0.15, rho=0.75){
  stopifnot(is.OptPx(o), o$Style$Quotient,is.OptPx(o), is.numeric(I1),is.numeric(I2),is.numeric(g1),is.numeric(g2),is.numeric(sigma1),
            is.numeric(sigma2), is.numeric(rho))
  
  sigma_a=sqrt(sigma1^2-2*rho*sigma1*sigma2+sigma2^2)
  d_ra12=(1/(sigma_a*sqrt(o$ttm)))*(log(I1/(o$K*I2))+(g2-g1-0.5*sigma1^2+0.5*sigma2^2)*o$ttm)
  d_1ra12=d_ra12+sigma_a*sqrt(o$ttm)
  G=(g2-g1-o$r)*o$ttm+sigma2*(sigma2-rho*sigma1)
  A=(I1/I2)*exp(G)
  B=stats::pnorm(d_1ra12)
  C=stats::pnorm(-d_1ra12)
  D=o$K*exp(-o$r*o$ttm)
  E=stats::pnorm(d_ra12)
  H=stats::pnorm(-d_ra12)
  
  # Calculate Price:
  
  if (o$Right$Name=='Call'){
    o$PxBS =A*B-D*E}
  
  if (o$Right$Name=='Put'){
    o$PxBS =-(A*C-D*H)}
    
  return (o)
}


#' @title Quotient option valuation via Monte Carlo (MC) model
#' @description Calculates the price of a Quotient option using Monte-Carlo simulations.
#' @details The Monte-Carlo simulations assume the underlying price undergoes Geometric Brownian Motion (GBM). 
#' Payoffs are discounted at risk-free rate to price the option.
#' A thorough understanding of the object class construction is recommended. 
#' Please see \code{OptPx}, \code{Opt} for more information. 
#' @author Richard Huang, Department of Statistics, Rice University, Spring 2015
#' 
#' @param o The \code{OptQuotient} option object to price. 
#' @param S0_2 The spot price of the second underlying asset.
#' @param NPaths Number of monte-carlo simulations to run. Larger number of trials lower variability at the expense of computation time.
#' 
#' @return An original \code{OptPx} object with Px.MC field as the price of the option and user-supplied S0_2, NPaths parameters attached.  
#'
#' @examples
#' (o = QuotientMC())$PxMC #Default Quotient option price. 
#' 
#' o = OptPx(Opt(S0=100, ttm=1, K=1.3), r=0.10, q=0, vol=0.1)
#' (o = QuotientMC(o, S0_2 = 180, NPaths=5))$PxMC
#' 
#' QuotientMC(OptPx(Opt()), S0_2 = 180, NPaths=5) 
#' 
#' QuotientMC(OptPx(), S0_2 = 201, NPaths = 5)
#' 
#' QuotientMC(OptPx(Opt(S0=500, ttm=1, K=2)), S0_2 = 1000, NPaths=5) 
#'  
#' @references 
#' \url{http://www.investment-and-finance.net/derivatives/q/quotient-option.html}
#'
#' @export
#' 
QuotientMC = function(o = OptPx(Opt(Style='Quotient')), S0_2= 100, NPaths = 5) {
  stopifnot(is.OptPx(o), is.numeric(NPaths), is.numeric(S0_2))  
  
  S0 = o$S0; vol = o$vol; r = o$r; ttm = o$ttm; Right = o$Right #Rename inputs for readability
  
  #Simulate stock price paths with GBM
  gbm <- function(S0=o$S0, vol=o$vol, r=o$r, ttm=o$ttm) 
    return(S0 * exp((r - 0.5 * vol^2) * ttm + vol * sqrt(ttm) * stats::rnorm(1)))
  
  #Binary payoff functions
  quot.payoff <- function(K= o$K, quot = 1) {
    if (o$Right$Call == TRUE) {
      return(max(quot-K, 0))
    }
    else if (o$Right$Put == TRUE){
      return(max(K - quot, 0))
    }
  }
  
  #Payoff vector
  payoff_vec = sapply(replicate(NPaths, gbm(S0 = S0_2)/gbm()), function(x) quot.payoff(quot=x))
  
  #Option price
  Px = exp(-r * ttm) * (sum(payoff_vec) / NPaths)
  
  #Return result
  o$S0_2 = S0_2
  o$NPaths = NPaths
  o$PxMC = Px  #add Binary MC option price
  return(o) 
}

