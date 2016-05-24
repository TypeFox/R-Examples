
#' @title DeferredPaymentLT
#' @description A binomial tree pricer of a Deferred Payment option. 
#' An American option that has payment at expiry no matter when exercise, 
#' causing differences in present value (PV) of a payoff.
#' @author Max Lee, Department of Statistics, Rice University, Spring 2015
#' 
#' @param o An object of class \code{OptPx}
#' @return An object of class \code{OptPx} with price included
#' 
#' @references Hull, J.C., \emph{Options, Futures and Other Derivatives}, 9ed, 2014. Prentice Hall. 
#' ISBN 978-0-13-345631-8, \url{http://www-2.rotman.utoronto.ca/~hull/ofod/index.html}
#' 
#' @examples
#' (o = DeferredPaymentLT())$PxLT
#' 
#' o = Opt(Style='DeferredPayment', Right="Call", S0=110,ttm=.5,K=110)
#' (o = DeferredPaymentLT(OptPx(o,r=.05,q=.04,vol=.2,NSteps=5)))$PxLT
#' 
#' o = Opt(Style='DeferredPayment', Right="Put", S0 = 50, ttm=2,K=47)
#' (o = DeferredPaymentLT(OptPx(o,r=.05,q=.04,vol=.25,NSteps=3)))$PxLT
#' 
#' @export
#' 
DeferredPaymentLT = function(o=OptPx(Opt(Style="DeferredPayment"))){
  stopifnot(o$Style$DeferredPayment, is.OptPx(o)); 
  
  u1 <- exp(o$vol*sqrt(o$ttm/o$NSteps))
  d1 <- 1/u1
  p1 <- (exp((o$r-o$q)*(o$ttm/o$NSteps))-d1)/(u1-d1)
  S = with(o, S0*d1^(0:NSteps)*u1^(NSteps:0)) 
  V = with(o, S0*d1^(0:NSteps)*u1^(NSteps:0))
  O = pmax(o$Right$SignCP * (S - o$K), 0) 
  
  
  ReCalc.O_S.on.Prior.Time.Step = function(i) {  
    O <<- if(i==o$NSteps){
      exp(-o$r * o$ttm/o$NSteps) * (p1*O[-i-1] + (1-p1)*O[-1])}else{
        exp(-o$r * o$ttm/o$NSteps*(o$NSteps-i)) * (p1*O[-i-1] + (1-p1)*O[-1])
      } 
    S <<- d1 * S[-i-1]                   
    if(i==o$NSteps){
      Payout = pmax(o$Right$SignCP *(S-o$K),0)
    }else{
      Payout = pmax(o$Right$SignCP *(S-o$K),0)*exp(-o$r*o$ttm/o$NSteps*(o$NSteps-i))
    }
    O <<- pmax(O, Payout)
    return(cbind(S, O))
  }
  
  BT = append(list(cbind(S, O)), sapply(o$NSteps:1, ReCalc.O_S.on.Prior.Time.Step)) 
  o$PxLT = BT[[length(BT)]][[2]] 
  return(o)
} 

