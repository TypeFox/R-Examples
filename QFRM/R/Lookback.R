#' @title Lookback option valuation with Black-Scholes (BS) model
#' @description Calculates the price of a lookback option using a BSM-adjusted algorithm; 
#' Carries the assumption that the asset price is observed continuously. 
#' @author Richard Huang, Department of Statistics, Rice University, Spring 2015
#' 
#' @param o An object of class \code{OptPx}.
#' @param Smax The maximum asset price observed to date.
#' @param Smin The minimum asset price observed to date.
#' @param Type Specifies the Lookback option as either Floating or Fixed- default argument is Floating.
#' 
#' @details To price the lookback option, we require the Smax/Smin, S0, r, q, vol, and ttm arguments 
#' from the object classes defined in the package. An example of a complete OptLookback option object can be found in the examples.
#' 
#' @return An original \code{OptPx} object with \code{PxBS} field as the price of the option 
#' and user-supplied \code{Smin}, \code{Smax}, and \code{Type} lookback parameters attached. 
#'
#' @references Hull, J.C., \emph{Options, Futures and Other Derivatives}, 9ed, 2014. Prentice Hall. 
#' ISBN 978-0-13-345631-8, \url{http://www-2.rotman.utoronto.ca/~hull/ofod/index.html}.
#' 
#' @examples
#'   (o = LookbackBS())$PxBS
#'   LookbackBS(OptPx(Opt(Style = 'Lookback'))) #Uses default arguments
#'   
#'   # See Hull 9e Example 26.2, p.608; gives price of 7.79 
#'   o = Opt(Style = 'Lookback', S0 = 50, ttm= .25, Right = "Put")
#'   o = OptPx(o,r = .1, vol = .4)
#'   o = LookbackBS(o, Type = "Floating") 
#'    
#'   # See Hull 9e Example 26.2, p.608; gives price of 8.04
#'   o = Opt(Style = 'Lookback', S0 = 50, ttm= .25, Right = "Call")
#'   o = OptPx(o, r = .1, vol = .4)
#'   o = LookbackBS(o, Type = "Floating") 
#'    
#'   # Price = 17.7129
#'   o = Opt(Style = 'Lookback', S0 = 50, ttm= 1, Right = "Put", K = 60)
#'   o = OptPx(o,r = .05, q = .02, vol = .25)
#'   o = LookbackBS(o, Type = "Fixed") 
#'
#'   # Price = 8.237
#'   o = Opt(Style = 'Lookback', S0 = 50, ttm= 1, Right = "Call", K = 55)
#'   o = OptPx(o,r = .1, q = .02, vol = .25)
#'   o = LookbackBS(o, Type = "Fixed")  
#'
#' @export
#' 
LookbackBS <- function(o=OptPx(Opt(Style='Lookback')), Smax = 50, Smin=50, Type = c("Floating", "Fixed")) {  
  stopifnot(is.OptPx(o), o$Style$Name=='Lookback', is.numeric(Smax), is.numeric(Smin), is.character(Type))  
  Type=match.arg(Type)
  isFixed = switch(Type, Fixed=TRUE, Floating=FALSE)
  o$Smax=Smax;o$Smin = Smin; o$Type = Type
  
  #Internal Pricing Functions
  fput.px <- function(S.max) {
    b1 = with(o, (log(S.max/S0) + (-r + q + vol^2/2)*ttm)/(vol*sqrt(ttm)))
    b2 = with(o, b1 - vol*sqrt(ttm))
    b3 = with(o, (log(S.max/S0) + (r - q - vol^2/2)*ttm)/(vol*sqrt(ttm)))
    Y2 = with(o, (2*(r - q - vol^2/2)*log(S.max/S0))/vol^2)  
    p_fl = with(o, S.max*exp(-r*ttm)*(pnorm(b1) - (vol^2*exp(Y2)*pnorm(-b3))/(2*(r-q))) + (S0*exp(-o$q*ttm)*vol^2*pnorm(-b2))/(2*(r-q)) - S0*exp(-q*ttm)*pnorm(b2))
    return(p_fl)
  }
  
  fcall.px <- function(S.min) {
    a1 = with(o,(log(S0/S.min) + (r - q + vol^2/2)*ttm)/(vol*sqrt(ttm)))
    a2 = with(o, a1 - vol*sqrt(ttm))
    a3 = with(o, (log(S0/S.min) + (-r + q + vol^2/2)*ttm)/(vol*sqrt(ttm)))
    Y1 =  with(o, (-2*(r - q - vol^2/2)*log(S0/S.min))/vol^2)  
    c_fl = with(o, S0*exp(-q*ttm)*pnorm(a1) - (S0*exp(-q*ttm)*vol^2*pnorm(-a1))/(2*(r-q)) - S.min*exp(-r*ttm)*(pnorm(a2) - (vol^2*exp(Y1)*pnorm(-a3))/(2*(r-q))))   
    return(c_fl)
  }
  
  #Logic Return Arguments
  if (o$Right$Name == "Put"){
    o$PxBS <- with(o, if (!isFixed) fput.px(Smin) else fcall.px(min(Smin, K)) + K*exp(-r*ttm) - S0*exp(-q*ttm))   
  } 
  else if (o$Right$Name == "Call") {
    o$PxBS <- with(o, if (!isFixed) fcall.px(Smax) else fput.px(max(Smax, K)) + S0*exp(-q*ttm) - K*exp(-r*ttm))  
  }
  
  return(o)
}


#' @title Lookback option valuation via Monte Carlo (MC) simulation
#' @description Calculates the price of a lookback option using a Monte Carlo (MC) Simulation.
#' Carries the assumption that the asset price is observed continuously.
#' Assumes that the the option o followes ds = mu * S * dt + sqrt(vol) * S * dz
#' where dz is a Wiener Process. Assume that without dividends, mu are default to be r. 
#'
#' @author Tong Liu, Department of Statistics, Rice University, Spring 2015
#' 
#' @param o The \code{OptPx} option object to price. See \code{OptPx} and \code{Opt} for more information.
#' @param NPaths How many time of the simulation are applied. Coustomer defined.
#' @param div number to decide length of each interval
#' @param Type Specifies the Lookback option as either Floating or Fixed- default argument is Floating.
#' 
#' @details To price the lookback option, we require the S0, K, and ttm arguments from object \code{Opt}
#' and r, q, vol from object OptPx defined in the package. The results of simulation would
#' unstable without setting seeds.
#' 
#' @return A list of class \code{LookbackMC} consisting of the input object \code{OptPx} and the price of the lookback option based on Monte Carlo Simulation (see references). 
#'
#' @examples
#'  (o = LookbackMC())$PxMC   #Use default arguments, Output: approximately 16.31.
#'   
#'  # Floating & Put
#'  o=OptPx(Opt(S0=50,K=50,ttm=0.25,Right='Put',Style="Lookback"),r=0.1,vol=.4)
#'  LookbackMC(o,NPaths=5,div=1000) #Output: 7.79 from Hull 9e Example 26.2 Pg 608. 
#'   
#'  # Floating & Call
#'  o=OptPx(Opt(S0=50,K=50,ttm=0.25,Right='Call',Style="Lookback"),r=0.1,vol=.4)
#'  LookbackMC(o,NPaths=5,div=1000) #Output: 8.04 from Hull 9e Example 26.2 Pg 608 
#'   
#'  # Fixed & Put
#'  o=OptPx(Opt(S0=50,K=60,ttm=1,Right='Put',Style="Lookback"),r=0.05,q=0.02,vol=.25)
#'  LookbackMC(o,Type="Fixed",NPaths=5,div=1000)
#'
#'  # Fixed & Call
#'  o=OptPx(Opt(S0=50,K=55,ttm=1,Right='Call',Style="Lookback"),r=0.1,vol=.25)
#'  LookbackMC(o,Type="Fixed",NPaths=5,div=1000)
#'  
#'  @references Hull, John C., \emph{Options, Futures and Other Derivatives}, 9ed, 2014. Prentice Hall. 
#'  ISBN 978-0-13-345631-8, \url{http://www-2.rotman.utoronto.ca/~hull/ofod}
#'
#' @export
#' 
LookbackMC=function(o=OptPx(Opt(Style='Lookback'), r=0.05, q=0, vol=0.3),
                    NPaths=5, div=1000, Type = c("Floating", "Fixed")){
  stopifnot(is.OptPx(o), o$Style$Name=='Lookback',is.character(Type), is.numeric(NPaths),is.numeric(div))  #---Condition
  o.class = class(o)
  #--- Is.Fixed and Is.Floating
  Type=match.arg(Type)
  isFixed = switch(Type, Fixed=T, Floating=F)
  isFloating = !isFixed
  
  dt=o$ttm/div      
  
  payoffs=replicate(NPaths,{
    ds.s=(o$r-o$q)*dt+o$vol*sqrt(dt)*stats::rnorm(div)
    st=c(o$S0,o$S0*cumprod(1+ds.s))
    
    smax=max(st)
    smin=min(st)
    
    if (isFloating) {
      K=smax*o$Right$Put+smin*o$Right$Call
      ST=st[length(st)]
    } else if (isFixed){
      K=o$K
      ST=smin*o$Right$Put+smax*o$Right$Call
    }
    payoff=exp(-o$r*o$ttm)*o$Right$SignCP*(max(ST-K))
  })
  p=mean(payoffs)
  
  o$isFixed=isFixed
  o$isFloating = !isFixed
  o$NPaths=NPaths
  o$PxMC=p
  
  class(o)=o.class   
  return(o)
}


