#' @title Average Strike option valuation via Monte Carlo (MC) simulation
#' @description Calculates the price of an Average Strike option using Monte Carlo simulations 
#' by determining the determine expected payout. Assumes that the input option follows a General
#' Brownian Motion \eqn{ds = mu * S * dt + sqrt(vol) * S * dz} where \eqn{dz ~ N(0,1)}
#' Note that the value of \eqn{mu} (the expected price increase) is assumped to be 
#' \code{o$r}, the risk free rate of return. Additionally, the averaging period is
#' assumed to be the life of the option. 
#' 
#'  @author Jake Kornblau, Department of Statistics and Department of Computer Science, Rice University, Spring 2015
#'  @param o The AverageStrike \code{OptPx} option to price. 
#'  @param NPaths the number of simulations to use in calculating the price,
#'  @return The original option object \code{o} with the price in the field \code{PxMC} based on the MC simulations. 
#'  
#'  @references Hull, John C., \emph{Options, Futures and Other Derivatives}, 9ed, 2014. Prentice Hall. 
#'  ISBN 978-0-13-345631-8, \url{http://www-2.rotman.utoronto.ca/~hull/ofod/index.html}
#'  Also, \url{http://www.math.umn.edu/~spirn/5076/Lecture16.pdf}
#'
#'  @examples
#'   (o = AverageStrikeMC())$PxMC   #Price =~ $3.6
#'   
#'   o = OptPx(o = Opt(Style='AverageStrike'), NSteps = 5)
#'   (o = AverageStrikeMC(o))$PxMC # Price =~ $6
#'   
#'   (o = AverageStrikeMC(NPaths = 20))$PxMC  #Price =~ $3.4
#'   
#'   o = OptPx(o = Opt(Style='AverageStrike'), NSteps = 5)
#'   (o = AverageStrikeMC(o, NPaths = 20))$PxMC  #Price =~ $5.6
#'   
#'   @export
AverageStrikeMC = function(o = OptPx(o=Opt(Style='AverageStrike')), NPaths = 5) {
  stopifnot(is.OptPx(o), is.numeric(NPaths), NPaths>0, o$Style$AverageStrike);
  
  o$PxMC = mean(
    sapply(
      (1:NPaths), 
      function(trial_num) {
        ds_div_S = with(o, exp((r - 0.5 * vol^2) * dt + vol * sqrt(dt) * rnorm(NSteps)))
        
        # ds is the product of a RV and the previous price. cumprod with S0 at
        # the beginning will accomplish this.
        prices = cumprod(c(o$S0, ds_div_S))
        prices = prices[2:(length(prices))]
        
        payoff = max(o$Right$SignCP * (prices[length(prices)] - mean(prices)), 0)
        return(exp(-o$r * o$ttm) * payoff)
      })
  )
  
  return(o)
}