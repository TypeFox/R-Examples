#' @title ForwardStart option valuation via Black-Scholes (BS) model
#' @description Compute the price of Forward Start options using BSM. 
#' A forward start option is a standard European option whose strike price is set equal to current asset price at some prespecified future date.
#' Employee incentive options are basically forward start option
#' @author Tongyue Luo, Department of Statistics, Rice University, Spring 2015
#' 
#' @param o an \code{OptPx} object including basic information of an option
#' @param tts Time to start of the option (in years)
#' @return The original user-supplied \code{OptPX} object 
#' with price field \code{PxBS} and any other provided user-supplied parameters.
#' 
#' @details A standard European option starts at a future time tts.
#' @examples
#' (o = ForwardStartBS())$PxBS
#' 
#' o = OptPx(Opt(Style='ForwardStart', Right='Put'))
#' (o = ForwardStartBS(o))$PxBS
#' 
#' @references Hull, John C., \emph{Options, Futures and Other Derivatives}, 9ed, 2014. Prentice Hall. 
#' ISBN 978-0-13-345631-8.\url{http://www-2.rotman.utoronto.ca/~hull/ofod/index.html}.
#' p.602
#' 
#' @export
#'
# ForwardStartBS = function(o = OptPx(Opt(Style='ForwardStart')), tts = 0.1){
#   
# }
 
ForwardStartBS = function(o = OptPx(Opt(Style='ForwardStart')), tts = 0.1){
  stopifnot(is.OptPx(o), o$Style$ForwardStart, is.numeric(tts), 0<=tts, tts<=o$ttm)
  
  Tdiff= o$ttm-tts
  mu = o$r - o$q
  q = o$q
  S0 = o$S0
  r = o$r
  #isCall = o$isCall
  
  d2= (mu-(o$vol^2)/2)*Tdiff/(o$vol*sqrt(Tdiff))
  d1= d2 + o$vol*sqrt(Tdiff)
  K=S0*exp(r*tts)
  
  C_N_d2 = stats::pnorm(d2, mean = 0, sd = 1)
  C_N_d1 = stats::pnorm(d1, mean = 0, sd = 1)
  c = S0*exp(-q*tts)*C_N_d1-K*exp(-r*Tdiff)*C_N_d2
  
  P_N_d1 = stats::pnorm((-d1), mean = 0, sd = 1)
  P_N_d2 = stats::pnorm((-d2), mean = 0, sd = 1)
  p = K * exp(-r * Tdiff) * P_N_d2 - S0 * exp(-q * tts) * P_N_d1
  
  o$tts = tts
  o$PxBS = (if (o$Right$Call) c else p) * exp(-q*tts)
  
  return(o)
}




#' @title Forward Start option valuation via Monte-Carlo (MC) simulation 
#' @description S3 object pricing model for a forward start European option using Monte Carlo simulation
#' @details A standard European option starts at a future time tts.
#' @author Tongyue Luo, Rice University, Spring 2015.
#' 
#' @param o An object of class \code{OptPx}
#' @param tts Time to start of the option, in years.
#' @param NPaths The number of MC simulation paths.
#' 
#' @return A list of class \code{ForwardStartMC} consisting of the input object
#'  \code{OptPx} and the appended new parameters and option price.
#' 
#' @references Hull, John C., \emph{Options, Futures and Other Derivatives}, 9ed, 2014. Prentice Hall. ISBN 978-0-13-345631-8, 
#' \url{http://www-2.rotman.utoronto.ca/~hull/ofod/index.html}.
#' \cr \url{http://investexcel.net/forward-start-options/}
#' 
#' @examples
#' (o = ForwardStartMC())$PxMC 
#' 
#' o = OptPx(Opt(Style='ForwardStart'), q = 0.03, r = 0.1, vol = 0.15)
#' (o = ForwardStartMC(o, tts=0.25))$PxMC 
#' 
#' ForwardStartMC(o = OptPx(Opt(Style='ForwardStart', Right='Put')))$PxMC
#' @export
#' 
ForwardStartMC = function (o= OptPx(Opt(Style='ForwardStart')), tts=0.1, NPaths = 5){
  stopifnot(is.OptPx(o), o$Style$ForwardStart, is.numeric(tts), 0<=tts, tts<=o$ttm, is.numeric(NPaths), 0<NPaths) 
  
#   T2 = o$ttm
#   S0 = o$S0
#   vol = o$vol
#   R = o$r
#   mu = o$r -o$q
#   Right = o$Right
  
#   isCall = switch(Right, "Call"=TRUE, "Put"=FALSE)
  
  CalcPrice = function (i){ 
    e = stats::rnorm(1)
    #Expectated stock price at tts (start time)
    S1 = with(o, S0 * exp (r*tts))
    S2 = with(o, S1 * exp((SYld - (vol^2)/2)*(ttm-tts) + vol* e* sqrt(ttm-tts)))
    
    if (o$Right$Call) {
      if (S2>S1) Payoff = S2 - S1
      else Payoff = 0
    }
    else {
      if (S2>S1) Payoff = 0
      else Payoff = S1 - S2
    }
    
    Price = exp(-o$r*o$ttm)*Payoff
    return(Price)
  }
  
  price = mean(sapply(1:NPaths, CalcPrice))
  o$PxMC = price
  return(o)
}



