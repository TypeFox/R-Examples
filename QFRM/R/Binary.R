#' @title Binary option valuation with Black-Scholes (BS) model
#' @description S3 object pricing model for a binary option.
#' Two types of binary options are priced: \code{'cash-or-nothing'} and \code{'asset-or-nothing'}.
#' @author Xinnan Lu, Department of Statistics, Rice University, Spring 2015
#'
#' @param o An object of class \code{OptPx}
#' @param Q A fixed amount of payoff
#' @param Type Binary option type: 'Cash or Nothing' or 'Asset or Nothing'. 
#' Partial names are allowed, eg. \code{'C'} or \code{'A'}
#' @return A list of class \code{Binary.BS} consisting of the input object \code{OptPx} and the appended new parameters and option price.
#'
#'  @references Hull, John C., \emph{Options, Futures and Other Derivatives}, 9ed, 2014. Prentice Hall. 
#' ISBN 978-0-13-345631-8, \url{http://www-2.rotman.utoronto.ca/~hull/ofod/index.html}.
#' pp.606-607
#' 
#' @examples
#' (o = BinaryBS())$PxBS
#' 
#' #This example should produce price 4.33 (see Derivagem, DG201.xls)
#' o = Opt(Style="Binary", Right='Call', S0=50, ttm=5/12, K=52)
#' o = OptPx(o, r=.1, vol=.40, NSteps=NA)
#' (o = BinaryBS(o, Q = 10, Type='cash-or-nothing'))$PxBS
#' 
#' BinaryBS(OptPx(Opt(Style="Binary"), q=.01), Type='asset-or-nothing')
#' BinaryBS(OptPx(Opt(Style="Binary", S0=100, K=80),q=.01))
#' o = Opt(Style="Binary", Right="Put", S0=50, K=60)
#' BinaryBS(OptPx(o,q=.04), Type='asset-or-nothing')
#' @export
#'
BinaryBS = function(o=OptPx(Opt(Style='Binary')), Q=1, Type=c('cash-or-nothing', 'asset-or-nothing')){
  stopifnot(is.OptPx(o), o$Style$Binary)
  
  Type=match.arg(Type)
  isCash = switch(Type, 'cash-or-nothing'=TRUE, 'asset-or-nothing'=FALSE)
  
  d1=BS(o)$BS$d1
  d2=BS(o)$BS$d2
  
  o = append(o, list(Q=Q, Type=Type, isCash=isCash, isAsset=!isCash)) 
  
  if (isCash) {
    c=with(o, Q*exp(-r*ttm)*pnorm(d2))
    p=with(o, Q*exp(-r*ttm)*pnorm(-d2))
  } else {
    c=with(o, S0*exp(-q*ttm)*pnorm(d1))
    p=with(o, S0*exp(-q*ttm)*pnorm(-d1))
  }
  
  if (o$Right$Call) o$PxBS = c
  if (o$Right$Put) o$PxBS = p
  
  return(o)
}



#' @title Binary option valuation via Monte-Carlo (via) simulation.
#' @author Tongyue Luo, Rice University, Spring 2015. 
#' @param o An \code{OptPx} object
#' @param Q A fixed numeric amount of payoff
#' @param Type Binary option type: \code{'cash-or-nothing'} or \code{'asset-or-nothing'}. 
#' @param NPaths The number of simulation paths to use in calculating the price
#' Partial names are allowed, eg. \code{'c'} or \code{'a'}
#' 
#' @details 
#' Two types of binary options are priced: \code{'cash-or-nothing'} and \code{'asset-or-nothing'}.
#' @references Hull, John C., \emph{Options, Futures and Other Derivatives}, 9ed, 2014. Prentice Hall. 
#' ISBN 978-0-13-345631-8, \url{http://www-2.rotman.utoronto.ca/~hull/ofod/index.html}.
#' pp.606-607. 
#' 
#' @return The original input object \code{o} with added parameters and option price \code{PxMC}
#' @examples
#' (o = BinaryMC())$PxMC
#' 
#' o = OptPx(Opt(Style="Binary"))
#' (o = BinaryMC(o, Type="cash"))$PxMC
#' 
#' o = OptPx(Opt(Style="Binary"),q=0.01)
#' (o = BinaryMC(o, Type="asset"))$PxMC
#' 
#' o = OptPx(Opt(Style="Binary", S0=100, K=80),q=0.01)
#' (o = BinaryMC(o, Type="cash"))$PxMC
#' 
#' o = OptPx(Opt(Style="Binary", Right="Put", S0=50, K=60),q=0.04)
#' (o = BinaryMC(o, Type="asset"))$PxMC
#' @export
BinaryMC = function(o=OptPx(Opt(Style='Binary')), Q=25, Type=c('cash-or-nothing', 'asset-or-nothing'), NPaths=5){
  stopifnot(is.OptPx(o), o$Style$Binary, is.numeric(Q), is.numeric(NPaths))
  
  Type = match.arg(Type)
  isCash = switch(Type, 'cash-or-nothing'=TRUE, 'asset-or-nothing'=FALSE)
  
  #This CalcLnSt is a function to generate random number and calculate the final price of the asset
  CalcLnSt = function (i=1) {
    e=stats::rnorm(1,mean=0, sd=1)
    deltalnS = (o$SYld - (o$vol^2)/2)*o$ttm + o$vol* e*sqrt(o$ttm)
    ST = o$S0 * exp(deltalnS)
    
    #Calculate the payoff 
    if (isCash) {
      if (o$Right$Call)
      {if (ST > o$K) Payoff = Q else Payoff = 0  } else
      {if (ST > o$K) Payoff = 0 else Payoff = Q }
    } 
    else {
      if (o$Right$Call) 
      {if (ST > o$K) Payoff = ST else Payoff = 0} else
      {if (ST > o$K) Payoff = 0  else Payoff = ST}
    }
    
    Price = Payoff * exp(-o$r *  o$ttm)
    return(Price)
  }
  
  o.Price = mean(sapply(1:NPaths, CalcLnSt))
  o$PxMC=o.Price 
  return (o)
}



#' @title Vectorized implementation of European option pricing algorithm
#' 
#' @param o \code{OptPx} object
#' @param Type Binary option type: \code{'cash-or-nothing'} or \code{'asset-or-nothing'}
#' @param Q A fixed amount of payoff 
#' 
#' @return \code{OptPx} object
#' 
#' @examples
#' (o = Binary_BOPM_Eu())$PxBT #basic example with defaultl values
#' 
#' (o = Binary_BOPM_Eu(Type='asset', Q=100))$PxBT
#' 
#' o = OptPx(o=Opt(Style='Binary'), r=0.05, q=0.02, rf=0.0, vol=0.30, NSteps=100)
#' (o = Binary_BOPM_Eu(o, Type='cash', Q=1000))$PxBT
#'
#' o = OptPx(o=Opt(Style='Binary'), r=0.15, q=0.01, rf=0.05, vol=0.35, NSteps=150)
#' (o = Binary_BOPM_Eu(o,Type='asset', Q=150))$PxBT
#' 
#' o = OptPx(o=Opt(Style='Binary'), r=0.025, q=0.001, rf=0.0, vol=0.10, NSteps=200)
#' (o = Binary_BOPM_Eu(o, Type='cash', Q=20))$PxBT
#' 
#' 
.Binary_BOPM_Eu = function(
  o=OptPx(Opt(Style='Binary')), Type=c('cash-or-nothing', 'asset-or-nothing'), Q=1000){
  stopifnot(is.OptPx(o), o$Style$Binary, is.numeric(Q))
  
  Type=match.arg(Type)
  isCash = switch(Type, 'cash-or-nothing'=TRUE, 'asset-or-nothing'=FALSE)
  isAsset = !isCash
  
  with(o, {
    S = S0*d^(NSteps:0)*u^(0:NSteps) # vector of terminal stock prices, lowest to highest (@t=ttm)
    
    if(isCash==1){
      O = vector('numeric',length=length(S))
      O = sapply(1:length(S), function(x) if(Right$SignCP*(S[x]-K) >= 0){O[x] = Q} else {O[x] = 0})
    } else if(isAsset==1){
      O = vector('numeric',length=length(S))
      O = sapply(1:length(S), function(x) if(Right$SignCP*(S[x]-K) >= 0){O[x] = S[x]} else {O[x] = 0})
    }
    
    csl = cumsum(log(c(1,1:NSteps)))       #-- logs avoid overflow & truncation
    tmp = csl[NSteps+1] - csl - csl[(NSteps+1):1] + log(p)*(0:NSteps) + log(1-p)*(NSteps:0)
    o$PxBT = DF_ttm * sum(exp(tmp)*O);     o$Type = Type; o$Q = Q
    return(o) #-- spot option price (at time 0)
  })
}



#' @title Binary option valuation vialattice tree (LT) implementation
#' @description Compute option price via binomial option pricing model (recombining symmetric binomial tree)
#' 
#' @param o \code{OptPx} object
#' @param Type Binary option type: \code{'cash-or-nothing'} or \code{'asset-or-nothing'}
#' @param Q A fixed amount of payoff
#' @param IncBT TRUE/FALSE, indicates whether to include the full binomial tree in the returned object
#' @return original \code{OptPx} object with \code{Px.BOPM} property and (optional) binomial tree
#'         IncBT = FALSE: option price value (type double, class numeric)
#'         IncBT = TRUE: binomial tree as a list (of length (o$n+1) of numeric matrices (2 x i). 
#'         Each matrix is a set of possible i outcomes at time step i 
#'         columns: (underlying prices, option prices)
#' 
#' @examples
#' (o = Binary_BOPM())$PxBT
#' 
#' o = OptPx(o=Opt(Style='Binary'))
#' (o = Binary_BOPM(o, Type='cash', Q=100, IncBT=TRUE))$PxBT
#' 
#' o = OptPx(Opt(Style='Binary'), r=0.05, q=0.02, rf=0.0, vol=0.30, NSteps=5)
#' (o = Binary_BOPM(o, Type='cash', Q=1000, IncBT=FALSE))$PxBT
#' 
#' o = OptPx(o=Opt(Style='Binary'), r=0.15, q=0.01, rf=0.05, vol=0.35, NSteps=5)
#' (o = Binary_BOPM(o,Type='asset',Q=150, IncBT=FALSE))$PxBT
#' 
#' o = OptPx(o=Opt(Style='Binary'), r=0.025, q=0.001, rf=0.0, vol=0.10, NSteps=5)
#' (o = Binary_BOPM(o, Type='cash', Q=20, IncBT=FALSE))$PxBT
#'
#' @export
#' 
Binary_BOPM = function(
  o=OptPx(Opt(Style='Binary')), Type=c('cash-or-nothing', 'asset-or-nothing'), Q=1000, IncBT=FALSE){ 
  
  stopifnot(is.OptPx(o), o$Style$Binary, is.numeric(Q)) # algorithm requires that a BinaryLT object is provided
  #NSteps=o$NSteps; p=o$p; K=o$K; S0 = o$S0
  
  Type=match.arg(Type);   isCash = ifelse(Type=='cash-or-nothing', TRUE, FALSE)
  
  if (!IncBT) return(.Binary_BOPM_Eu(o)) else { 
    S = with(o, o$S0*d^(0:o$NSteps)*u^(o$NSteps:0)) # vector of terminal stock prices, lowest to highest (@t=ttm)
    
    if(isCash){
      O = vector('numeric',length=length(S))
      O = sapply(1:length(S), function(x) if(o$Right$SignCP*(S[x]-o$K) >= 0){O[x] = Q} else {O[x] = 0})
    } else {
      O = vector('numeric',length=length(S))
      O = sapply(1:length(S), function(x) if(o$Right$SignCP*(S[x]-o$K) >= 0){O[x] = S[x]} else {O[x] = 0})
    }
    
    RecalcOSonPriorTimeStep = function(i) { #sapply(1:(i-1), function(j) 
      O <<- o$DF_dt * (o$p*O[-i-1] + (1-o$p)*O[-1])  #prior option prices (@time step=i-1)
      S <<- o$d * S[-i-1]                   # prior stock prices (@time step=i-1)
      Payout = pmax(o$sgnCP * (S - o$K), 0)   # payout at time step i-1 (moving backward in time)
      if (o$Style$American) O <<- pmax(O, Payout)    # 
      return(cbind(S, O))
    }
    
    BT = append(list(cbind(S, O)), sapply(o$NSteps:1, RecalcOSonPriorTimeStep)) #binomial tree
    o$PxBT = BT[[length(BT)]][[2]]; o$Type=Type; o$Q = Q # add BOPM price
    if (IncBT) o$BT = BT
    return(o)
  } 
}
