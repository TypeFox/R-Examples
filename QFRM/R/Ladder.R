#' @title Ladder option valuation via Monte Carlo (MC) simulation.
#' @description Calculates the price of a Ladder Option using 5000 Monte Carlo simulations.
#' The helper function LadderCal() aims to calculate expected payout for each stock prices.
#' 
#' \emph{Important Assumptions}:
#' The option o follows a General Brownian Motion (BM) 
#' \eqn{ds = mu * S * dt + sqrt(vol) * S * dW} where \eqn{dW ~ N(0,1)}. 
#' The value of \eqn{mu} (the expected price increase) is assumed to be \code{o$r}, the risk free rate of return. 
#' 
#' @author Huang Jiayao, Risk Management and Business Intelligence at Hong Kong University of Science and Technology, 
#'          Exchange student at Rice University, Spring 2015
#'          
#' @param o The \code{OptPx} Ladder option object to price. 
#' @param NPaths The number of simulation paths to use in calculating the price
#' @param L A series of ladder strike price.
#' @return The option \code{o} with the price in the field \code{PxMC} based on MC simulations 
#'    and the ladder strike price \code{L} set by the users themselves 
#'  
#' @references  
#' \url{http://stackoverflow.com/questions/25946852/r-monte-carlo-simulation-price-path-converging-volatility-issue}
#'
#' @examples
#'  (o = LadderMC())$PxMC #Price = ~12.30
#'   
#'  o = OptPx(o=Opt(Style='Ladder'), NSteps = 5)
#'  (o = LadderMC(o))$PxMC        #Price = ~11.50
#'   
#'  o = OptPx(Opt(Style='Ladder', Right='Put'))
#'  (o = LadderMC(o, NPaths = 5))$PxMC   # Price = ~12.36
#'   
#'  (o = LadderMC(L=c(55,65,75)))$PxMC   # Price = ~10.25
#' @export
LadderMC = function(o = OptPx(o=Opt(Style='Ladder'), NSteps=5), NPaths = 5, L=c(60,80,100)) {
  stopifnot(is.OptPx(o), is.numeric(NPaths), NPaths>0, o$Style$Ladder);
  
  LadderCal = function(price)
  {
    
    #Purpose: The function aims to calculate payout for ladder options
    
    K=o$K
    i=length(L)
    
    downladder= max(price-K,0)
    
    upladder = max(sapply(L,function(n){if(price>n) n=n-K else n=0}))
    
    if (upladder == 0) payout = downladder else payout = upladder  
    return(payout)
  }
  
  o$L=L
  
  o$PxMC = o$Right$SignCP*mean(sapply((1:NPaths),function(trial_sum){
    div = with(o, exp((r - 0.5 * vol^2) * dt + vol * sqrt(dt) * rnorm(NSteps)))
    inprices = cumprod(c(o$S0, div))
    highestprice = max(inprices)
    Ladderprice = LadderCal(highestprice)
    return(exp(-o$r * o$ttm) * Ladderprice)
  }))
  
  return(o)
}
