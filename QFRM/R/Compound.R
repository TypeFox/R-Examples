#' @title Compound option valuation with Black-Scholes (BS) model
#' @author Robert Abramov
#' 
#' @param o = \code{OptPx} object
#' @param K1 The first Strike Price (of the option on the option)
#' @param T1 The time of first expiry (of the option on the option)
#' @param Type Possible choices are 
#' \code{cc} - call option on call option
#' \code{cp} - call on put 
#' \code{pc} - put on call 
#' \code{pp} - put on put 
#' 
#' @return A list of object 'OptCompound' containing the option parameters binomial tree parameters and compound option parameters
#' 
#' @examples
#' (o <- CompoundBS())$PxBS #price compound option with default parameters 
#' 
#' o = OptPx(Opt(Style='Compound'), r=0.05, q=0.0, vol=0.25)
#' CompoundBS(o,K1=10,T1=0.5)
#' 
#' o = Opt(Style='Compound', S0=50, K=52, ttm=1)
#' CompoundBS(o=OptPx(o, r=.05, q=0, vol=.25),K1=6,T1=1.5)
#' 
#' o = Opt(Style='Compound', S0=90, K=100, ttm=1.5)
#' CompoundBS(o=OptPx(o, r=.05, q=0, vol=.25),K1=15,T1=1)
#' 
#' o = Opt(Style='Compound', S0=15, K=15, ttm=0.25)
#' CompoundBS(o=OptPx(o, r=.05, q=0, vol=.25),K1=3,T1=1.5)
#' 
#' @export
CompoundBS = function(o=OptPx(Opt(Style='Compound')), K1 = 10, T1 = 0.5, Type=c('cc','cp','pp','pc')){
  stopifnot(o$Style$Compound, is.OptPx(o), is.numeric(K1), is.numeric(T1))
  K2 = o$K
  S0 = o$S0  
  T2 = T1 + o$ttm
  r = o$r
  q = o$q
  vol = o$vol
  
  Type = match.arg(Type)
  Sstar = stats::uniroot(function(x) BS_Simple(x,K2,r,q,T2-T1,vol)$Px$Call-K1, interval=c(0,1000))$root
  
  a1 = (log(S0/Sstar) + (r - q + .5*(vol^2))*T1)/(vol*sqrt(T1))
  a2 = a1 - vol*sqrt(T1)
  
  b1 = (log(S0/K2) + (r - q + .5*(vol^2))*T2)/(vol*sqrt(T2))
  b2 = b1 - vol*sqrt(T2)
  
  cc = max((S0*exp(-q*T2)*pbnorm(a1,b1,sqrt(T1/T2))) - (K2*exp(-r*T2)*pbnorm(a2,b2,sqrt(T1/T2))) - (exp(-r*T1)*K1*stats::pnorm(a2)),0)
  cp = max((K2*exp(-r*T2)*pbnorm(-a2,-b2,sqrt(T1/T2))) - (S0*exp(-q*T2)*pbnorm(-a1,-b1,sqrt(T1/T2))) - (exp(-r*T1)*K1*stats::pnorm(-a2)),0)
  pc = max((K2*exp(-r*T2)*pbnorm(-a2,b2,-sqrt(T1/T2))) - (S0*exp(-q*T2)*pbnorm(-a1,b1,-sqrt(T1/T2))) + (exp(-r*T1)*K1*stats::pnorm(-a2)),0)
  pp = max((S0*exp(-q*T2)*pbnorm(a1,-b1,-sqrt(T1/T2))) - (K2*exp(-r*T2)*pbnorm(a2,-b2,-sqrt(T1/T2))) + (exp(-r*T1)*K1*stats::pnorm(a2)),0)
  
  #o.class = class(o)   # remember the class name
  Px = list(CallonCall=cc, CallonPut=cp, PutonCall=pc, PutonPut=pp)
  BS = list(Sstar=Sstar, a1=a1, a2=a2, b1=b1, b2=b2, Px=Px)
  
  o$BS = BS
  o$K1 = K1
  o$T1 = T1
  o$Type = Type
  o$PxBS = if (Type=='cc') cc else {if (Type=='cp') cp else {if (Type=='pc') pc else {if (Type=='pp') pp else NA}}}
  
  #class(o)=c(o.class, 'OptCompound')   #-- Child 'Compound' inherits properties of parent 'OptBT'
  return(o)  
}



#' @title Compound option valuation via lattice tree (LT) model
#' @description \code{CompoundLT} prices a compound option using the binomial tree (BT) method. 
#' The inputs it takes are two \code{OptPx} objects. 
#' It pulls the S from the o2 input which should be the option with the greater time to maturity.
#' 
#' @author Kiryl Novikau, Department of Statistics, Rice University, Spring 2015
#' @param o1 The \code{OptPx} object with the shorter time to maturity 
#' @param o2 The \code{OptPx} object with the longer time to maturity
#'        
#' @return User-supplied \code{o1} option with fields \code{o2} and \code{PxLT}, 
#' as the second option and calculated price, respectively.
#' 
#' @references Hull, John C., \emph{Options, Futures and Other Derivatives}, 9ed, 2014. Prentice Hall. 
#' ISBN 978-0-13-345631-8, \url{http://www-2.rotman.utoronto.ca/~hull/ofod/index.html}.
#' 
#' @examples
#' (o = CompoundLT())$PxLT # Uses default arguments
#' 
#' #Put option on a Call:
#' o = Opt(Style="Compound", S0=50, ttm=.5, Right="P", K = 50)
#' o1 = OptPx(o, r = .1, vol = .4, NSteps = 5)
#' o = Opt(Style="Compound", S0=50, ttm=.75, Right="C", K = 60)
#' o2 = OptPx(o, r = .1, vol = .4, NSteps = 5)
#' (o = CompoundLT(o1, o2))$PxLT
#' 
#' #Call option on a Call:
#' o = Opt(Style = "Compound", S0 = 50, ttm= .5, Right = "Call", K = 50)
#' o1 = OptPx(o, r = .1, vol = .4, NSteps = 5)
#' o = Opt(Style = "Compound", S0 = 50, ttm= .75, Right = "Call", K = 5)
#' o2 = OptPx(o, r = .1, vol = .4, NSteps = 5)
#' (o = CompoundLT(o1, o2))$PxLT
#' 
#' #Put option on a Put:
#' o = Opt(Style = "Compound", S0 = 50, ttm= .5, Right = "Put", K = 40)
#' o1 = OptPx(o, r = .1, vol = .4, NSteps = 5)
#' o = Opt(Style = "Compound", S0 = 50, ttm= .75, Right = "Put", K = 50)
#' o2 = OptPx(o, r = .1, vol = .4, NSteps = 5)
#' (o = CompoundLT(o1, o2))$PxMC
#' 
#' #Call option on a Put:
#' o = Opt(Style = "Compound", S0 = 50, ttm= .5, Right = "Call", K = 30)
#' o1 = OptPx(o, r = .1, vol = .4, NSteps = 5)
#' o = Opt(Style = "Compound", S0 = 50, ttm= .75, Right = "Put", K = 80)
#' o2 = OptPx(o, r = .1, vol = .4, NSteps = 5)
#' (o = CompoundLT(o1, o2))$PxLT
#' 
#' @export
#' 
CompoundLT <- function(o1 = OptPx(Opt(Style = "Compound")), o2 = OptPx(Opt(Style = "Compound"))){
  
  stopifnot(o1$Style$Compound, o2$Style$Compound, is.OptPx(o1),is.OptPx(o1))
  #' Getting the necessary variables
  options(expressions = 5e5)
  K1 = o1$K; t1 = o1$ttm; K2 = o2$K; t2 = o2$ttm; z = o1$vol; r = o1$r
  q = o1$q; S = o2$S0; NSteps = o1$NSteps; dt = t2/ NSteps; u = exp(z * dt^(.5))
  d = 1 / u; a = exp((r - q)*dt); p = (a-d)/(u-d)
  sCP1 = o1$Right$Call*2-1; sCP2 = o2$Right$Call*2-1; temp = round(t1/dt)
  payouts <- pmax(sCP2*(S*d^(NSteps:0)*u^(0:NSteps) - K2),0)
  payouts <- payouts[length(payouts):1]
  payouts2 <- NA
  compound <- NA
  
  compoundfunc <- function(inp = (p*payouts[1:(NSteps+1)] + (1-p)*payouts[1:(NSteps+1)])){
    if(length(inp) == temp + 1){
      payouts2 <<- inp
    } else {
      Recall(inp = (inp[-length(inp)]*(p) + inp[-1]*(1-p))*exp(-r*dt))
    }
  }
  compoundfunc()
  
  #' Pricing
  compound <- pmax(sCP1*(payouts2 - K1), 0)
  csl = cumsum(log(c(1,1:temp)))   
  tmp = csl[temp+1] - csl - csl[(temp+1):1] + log(1-p)*(0:temp) + log(p)*(temp:0)
  pricefinal = exp(-r*t1) * sum(exp(tmp)*compound)
  o1$PxLT = pricefinal
  o1$o2 = o2
  return(o1)
}