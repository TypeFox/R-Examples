#' @title Asian option valuation via Black-Scholes (BS) model
#' @description Price Asian option using BS model 
#' @author Xinnan Lu, Department of Statistics, Rice University, Spring 2015
#' @details This pricing algorithm assumes average price is calculated continuously.
#' 
#' @param o An object of class \code{OptPx}
#' 
#' @return A list of class \code{AsianBS} consisting of the original \code{OptPx} object
#' and the option pricing parameters \code{M1}, \code{M2}, \code{F0}, and \code{sigma} 
#' as well as the computed option price \code{PxBS}.
#' 
#' @references Hull, John C., \emph{Options, Futures and Other Derivatives}, 9ed, 2014. Prentice Hall. 
#' ISBN 978-0-13-345631-8, \url{http://www-2.rotman.utoronto.ca/~hull/ofod/index.html} pp.609-611.
#' 
#' @examples
#'  (o = AsianBS())$PxBS #Price = ~4.973973,  using default values
#'  
#'  o = Opt(Style='Asian',S0=100,K=90,ttm=3)
#'  (o = AsianBS(OptPx(o,r=0.03,q=0,vol=0.3)))$PxBS
#'  
#'  o = Opt(Style='Asian',Right='P',S0=100,K=110,ttm=0.5)
#'  (o = AsianBS(OptPx(o,r=0.03,q=0.01,vol=0.3)))$PxBS
#'  
#'  #See J.C.Hull, OFOD'2014, 9-ed, ex.26.3, pp.610. The price is 5.62.
#'  o = Opt(Style='Asian',Right='Call',S0=50,K=50,ttm=1)
#'  (o = AsianBS(OptPx(o,r=0.1,q=0,vol=0.4)))$PxBS
#' @export
AsianBS = function( o=OptPx(Opt(Style='Asian'))){
  stopifnot(is.OptPx(o), o$Style$Asian)
  
  # calcualte the two moments M1, M2 of average price
  M1=with(o, (exp((r-q)*ttm)-1)*S0/((r-q)*ttm))
  M2=with(o, 2*S0^2*exp((2*(r-q)+vol^2)*ttm)/((r-q+vol^2)*(2*r-2*q+vol^2)*ttm^2)+
            2*S0^2/((r-q)*ttm^2)*(1/(2*(r-q)+vol^2)-exp((r-q)*ttm)/(r-q+vol^2)))
  
  F0=M1 # forward price
  sigma=with(o, sqrt((1/ttm)*log(M2/M1^2))) # implied volatility
  
  d1=with(o, (log(F0/K)+sigma^2*ttm/2)/(sigma*sqrt(ttm)))
  d2=with(o, d1-sigma*sqrt(ttm))
  
  o$M1=M1; o$M2=M2; o$F0=F0; o$sigma=sigma;
  o$PxBS=with(o, exp(-r*ttm)*(Right$SignCP*F0*pnorm(Right$SignCP*d1)-Right$SignCP*K*pnorm(Right$SignCP*d2)))
  
  return(o)
}




#' @title Asian option valuation with Monte Carlo (MC) simulation.
#' @description Calculates the price of an Asian option using Monte Carlo simulations to 
#' determine expected payout. 
#' \cr Assumptions:
#' \cr The option follows a General Brownian Motion (BM), 
#' \cr \eqn{ds = mu * S * dt + sqrt(vol) * S * dW} where \eqn{dW ~ N(0,1)}.
#' \cr The value of \eqn{mu} (the expected price increase) is \code{o$r}, the risk free rate of return (RoR). 
#' \cr The averaging period is the life of the option. 
#' 
#'  @author Jake Kornblau, Department of Statistics and Department of Computer Science, Rice University, 2016
#'  @param o The \code{OptPx} Asian option to price. 
#'  @param NPaths The number of simulation paths to use in calculating the price,
#'  @return The option \code{o} with the price in the field \code{PxMC} based on MC simulations. 
#'  
#'  @references Hull, John C., \emph{Options, Futures and Other Derivatives}, 9ed, 2014. 
#'    Prentice Hall. ISBN 978-0-13-345631-8, 
#'     \cr \url{http://www-2.rotman.utoronto.ca/~hull/ofod/index.html}
#'     \cr \url{http://www.math.umn.edu/~spirn/5076/Lecture16.pdf}
#'     
#'  @examples
#'   (o = AsianMC())$PxMC #Price = ~5.00,  using default values
#'   
#'   o = OptPx(Opt(Style='Asian'), NSteps = 5)
#'   (o = AsianMC(o, NPaths=5))$PxMC #Price = ~$5
#'   
#'   (o = AsianMC(NPaths = 5))$PxMC # Price = ~$5
#'   
#'   o = Opt(Style='Asian', Right='Put',S0=10, K=15)
#'   o = OptPx(o, r=.05, vol=.1, NSteps = 5)
#'   (o = AsianMC(o, NPaths = 5))$PxMC # Price = ~$4
#'   
#'   #See J.C.Hull, OFOD'2014, 9-ed, ex.26.3, pp.610. 
#'  o = Opt(Style='Asian',S0=50,K=50,ttm=1)
#'  o = OptPx(o,r=0.1,q=0,vol=0.4,NSteps=5)
#'  (o = AsianBS(o))$PxBS   #Price is 5.62.
#'  (o = AsianMC(o))$PxMC
#'   @export
#'   
AsianMC = function(o = OptPx(o=Opt(Style='Asian'), NSteps=5), NPaths = 5) {
  stopifnot(is.OptPx(o), is.numeric(NPaths), NPaths>0, o$Style$Asian);
  
  o$PxMC = mean(
    sapply(
      (1:NPaths), 
      function(trial_num) {
        ds_div_S = with(o, exp((r - 0.5 * vol^2) * dt + vol * sqrt(dt) * rnorm(NSteps)))
        
        # ds is the product of a RV and the previous price. cumprod with S0 at
        # the beginning will accomplish this.
        prices = cumprod(c(o$S0, ds_div_S))
        prices = prices[2:(length(prices))]
        
        payoff = max(o$Right$SignCP * (mean(prices) - o$K), 0)
        
        return(exp(-o$r * o$ttm) * payoff)
      })
  )
  
  return (o)
}