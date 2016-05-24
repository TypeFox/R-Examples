#' @title Gap option valuation via Black-Scholes (BS) model
#' @description S3 object constructor for price of gap option using BS model
#' @author Tong Liu, Department of Statistics, Rice University, Spring 2015
#' 
#' @param o An object of class \code{OptPx}
#' @param K2 Strike price that determine if the option pays off.
#' 
#' @return An original \code{OptPx} object with \code{PxBS} field as the price of the option 
#' and user-supplied \code{K2} parameter  
#' 
#' @references Hull, John C., \emph{Options, Futures and Other Derivatives}, 9ed, 2014. Prentice Hall. 
#' ISBN 978-0-13-345631-8. \url{http://www.mathworks.com/help/fininst/gapbybls.html}
#' 
#' @examples
#' #See J.C.Hull, OFOD'2014, 9-ed, Example 26.1, p.601
#' (o <- GapBS())$PxBS
#' 
#' GapBS(o=OptPx(Opt(Style='Gap',Right='Put',K=57)))    
#' 
#' #See http://www.mathworks.com/help/fininst/gapbybls.html
#' o = Opt(Style='Gap',Right='Put',K=57,ttm=0.5,S0=52)
#' o = GapBS(OptPx(o,vol=0.2,r=0.09),K2=50)
#' 
#' o = Opt(Style='Gap',Right='Put',K=57,ttm=0.5,S0=50)
#' (o <- GapBS(OptPx(o,vol=0.2,r=0.09),K2=50))$PxBS
#' @export
#' 
GapBS = function(
  o = OptPx(Opt(Style='Gap', Right='Put', S0=500000, K=400000, ttm=1, 
                ContrSize=1, SName='Insurance coverage example #26.1, p.601, OFOD, J.C.Hull, 9ed.'),
            r=.05,q=0,vol=.2), 
  K2 = 350000){
#   GapBS=function(o=OptPx(Opt(Style='Gap',K=57),vol=0.2,r=0.05,q=0),K2=50){
  stopifnot(is.OptPx(o), o$Style$Gap)
  
  d1=with(o, (log(S0/K2)+(r-q+vol^2/2)*ttm)/(vol*sqrt(ttm)))
  d2=d1-o$vol*sqrt(o$ttm)
  
  o$PxBS= with(o, Right$SignCP*(S0*exp(-q*ttm)*pnorm(Right$SignCP*d1)-K*exp(-r*ttm)*pnorm(Right$SignCP*d2)))
  o$K2=K2
  return(o)
}



#' @title Gap option valuation via lattice tree (LT) model
#' @description A binomial tree pricer of Gap options that takes the average results for given step sizes in NSteps.  
#' Large step sizes should be used for optimal accuracy but may take a minute or so.
#' 
#' @author Max Lee, Department of Statistics, Rice University, Spring 2015
#' 
#' @param o An object of class \code{OptPx}
#' @param K2 A numeric strike price above used in calculating if option is in the money or not, known as trigger.
#' @param on A vector of number of steps to be used in binomial tree averaging, vector of positive intergers.
#' @return An onject of class \code{OptPx} including price
#' 
#' @references Hull, John C., \emph{Options, Futures and Other Derivatives}, 9ed, 2014. Prentice Hall. 
#' ISBN 978-0-13-345631-8. \url{http://www-2.rotman.utoronto.ca/~hull/ofod/index.html}.
#' \cr Humphreys, Natalia. University of Dallas. 
#' 
#' @examples
#' (o = GapLT())$PxLT
#' 
#' o = Opt(Style="Gap",Right='Put',S0 = 500000, ttm = 1,K = 400000)
#' o = OptPx(o,r = .05, q=0, vol =.2)
#' (o = GapLT(o,K2 = 350000,on=c(498,499,500,501,502)))$PxLT
#'
#' o = Opt(Style="Gap", Right='Call',S0 = 65, ttm = 1,K = 70)
#' o = OptPx(o,r = .05, q=.02,vol =.1)
#  (o = GapLT(o,K2 = 60, on=c(498,499,500,501,502)))$PxLT
#' 
#' @export
#' 
GapLT = function(o=OptPx(Opt(Style="Gap")), K2=60, on=c(100,200)){
  
  stopifnot(o$Style$Gap, is.OptPx(o), is.numeric(K2), K2>0,
            is.vector(on), length(on)>=2, on[1]>0, on[2]>0)
  Px = c()
  
  Px = sapply(1:length(o$NSteps), function(x){
    u1 = exp(o$vol*sqrt(o$ttm/on[x]))
    d1 = 1/u1
    p1 = (exp((o$r-o$q)*(o$ttm/on[x]))-d1)/(u1-d1)
    
    with(o, {
      S = S0*d1^(on[x]:0)*u1^(0:on[x]) 
      
      O = vector('numeric',length=length(S))
      O = sapply(1:length(S), function(y) 
        if(o$Right$SignCP*(S[y]-K2) > 0) 
          {O[y] = Right$SignCP*(S[y]-o$K)} else {O[y] = 0})
      
      csl = cumsum(log(c(1,1:on[x])))       
      tmp = csl[on[x]+1] - csl - csl[(on[x]+1):1] + log(p1)*(0:on[x]) + log(1-p1)*(on[x]:0)
      Px[x] = exp(r*-ttm) * sum(exp(tmp)*O)
    })
  }) 
  
  o$PxLT = mean(Px)
  o$K2 = K2
  o$on = on
  return(o)
}



#' @title Gap option valuation via Monte Carlo (MC) simulation
#' @description GapMC prices a gap option using the MC method. 
#'   The call payoff is \eqn{S_T-K} when \eqn{S_T>K2}, where \eqn{K_2} is the trigger strike. 
#'   The payoff is increased by \eqn{K_2-K}, which can be positive or negative.
#'   The put payoff is \eqn{K-S_T} when \eqn{S_T<K_2}. 
#'   Default values are from policyholder-insurance example 26.1, p.601, from referenced OFOD, 9ed, text. 
#' @author Kiryl Novikau, Department of Statistics, Rice University, Spring 2015
#' 
#' @param o The \code{OptPx} object (See \code{OptPx()} constructor for more information)
#' @param K2 The trigger strike price.
#' @param NPaths The number of paths (trials) to simulate.
#' @return An \code{OptPx} object. The price is stored under \code{o$PxMC}.
#' 
#' @references Hull, John C., \emph{Options, Futures and Other Derivatives}, 9ed, 2014. Prentice Hall. 
#' ISBN 978-0-13-345631-8. \url{http://www-2.rotman.utoronto.ca/~hull/ofod/index.html}.
#' p.601
#' 
#' @examples
#' (o = GapMC())$PxMC  #example 26.1, p.601
#' 
#' o = Opt(Style='Gap', Right='Call', S0=50, K=40, ttm=1)
#' o = OptPx(o, vol=.2, r=.05, q = .02)
#' (o = GapMC(o, K2 = 45, NPaths = 5))$PxMC
#'
#' o = Opt(Style='Gap', Right='Call', S0 = 50, K = 60, ttm = 1)
#' o = OptPx(o, vol=.25,r=.15, q = .02)
#' (o = GapMC(o, K2 = 55, NPaths = 5))$PxMC
#' 
#' o = Opt(Style='Gap', Right = 'Put', S0 = 50, K = 57, ttm = .5)
#' o = OptPx(o, vol = .2, r = .09, q = .2)
#' (o = GapMC(o, K2 = 50, NPaths = 5))$PxMC
#' 
#' o = Opt(Style='Gap', Right='Call', S0=500000, K=400000, ttm=1)
#' o = OptPx(o, vol=.2,r=.05, q = 0)
#' (o = GapMC(o, K2 = 350000, NPaths = 5))$PxMC
#' 
#' @export
#' 
GapMC = function(
  o = OptPx(Opt(Style='Gap', Right='Put', S0=500000, K=400000, ttm=1, 
                ContrSize=1, SName='Insurance coverage example #26.1, p.601, OFOD, J.C.Hull, 9ed.'),
                           r=.05,q=0,vol=.2), 
                 K2 = 350000, NPaths = 5){
  stopifnot(o$Style$Gap, is.OptPx(o), is.numeric(K2), is.numeric(NPaths))
  
  gbm = function(S0 = o$S0, vol = o$vol, r = o$r, ttm = o$ttm, q = o$q){
    return(S0 * exp((r - q - (.5 * vol^2)) * ttm + vol * sqrt(ttm) * stats::rnorm(1)))
  }
  
  PayoffGap = function(S) (pmax(o$Right$SignCP*(S - K2)) > 0) * (o$Right$SignCP*(S - o$K))
  
  payoff = sapply(replicate(NPaths, gbm()), PayoffGap)
  Px = exp(-o$r * o$ttm) * (sum(payoff) / NPaths)

  o$K2 = K2
  o$NPaths = NPaths
  o$PxMC = Px
  return(o)
}

