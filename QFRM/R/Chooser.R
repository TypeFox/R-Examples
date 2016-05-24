#'@title       Chooser option valuation via Black-Scholes (BS) model 
#'@description Compute an exotic option that allow the holder decide the option 
#'              will be a call or put option at some predetermined future date.
#'              In a simple case, both put and call option are plain vanilla option.
#'      The value of the simple chooser option is \eqn{\max{C(S,K,t_1),P(S,K,t_2)}}.
#'              The plain vanilla option is calculated based on the BS model.
#'@author      Le You, Department of Statistics, Rice University, spring 2015
#'              
#'@param        o An object of class \code{OptPx}
#'@param        t1    The time to maturity of the call option, measured in years.
#'@param        t2    The time to maturity of the put option, measured in years.
#'@return       A list of class \code{SimpleChooserBS} consisting of the original \code{OptPx} object 
#'              and the option pricing parameters \code{t1}, \code{t2}, 
#'              as well as the computed price \code{PxBS}.
#'              
#'@references   
#'\itemize{
#' \item{Hull, John C.,\emph{Options, Futures and Other Derivatives}, 9ed, 2014. Prentice Hall. 
#'      ISBN 978-0-13-345631-8. \url{http://www-2.rotman.utoronto.ca/~hull/ofod/index.html}}
#' \item{Huang Espen G., \emph{Option Pricing Formulas}, 2ed. 
#'     \url{http://down.cenet.org.cn/upfile/10/20083212958160.pdf} }
#' \item{Wee, Lim Tiong, MFE5010 \emph{Exotic Options,Notes for Lecture 4 Chooser option}.
#'     \url{http://www.stat.nus.edu.sg/~stalimtw/MFE5010/PDF/L4chooser.pdf} }
#' \item{Humphreys, Natalia A., ACTS 4302 Principles of Actuarial Models: Financial Economics.
#'     \emph{Lesson 14: All-or-nothing, Gap, Exchange and Chooser Options}. 
#'      } }
#'              
#'@examples 
#'(o = ChooserBS())$PxBS
#'
#'o = Opt(Style='Chooser',Right='Other',S0=50, K=50)
#'(o = ChooserBS(OptPx(o, r=0.06, q=0.02, vol=0.2),9/12, 3/12))$PxBS
#'  
#'o = Opt(Style='Chooser',Right='Other',S0=50, K=50)
#'(o = ChooserBS (OptPx(o,r=0.08, q=0, vol=0.25),1/2, 1/4))$PxBS  
#'
#'o = Opt(Style='Chooser',Right='Other',S0=100, K=50)
#'(o = ChooserBS(OptPx(o,r=0.08, q=0.05, vol=0.3),1/2, 1/4))$PxBS
#'
#'@export 
#'
ChooserBS =function(o=OptPx(Opt(Style='Chooser')),t1=9/12,t2=3/12)
{ stopifnot(is.OptPx(o), o$Style$Chooser)
  d2 = (log(o$S0/o$K) + ((o$r-o$q) - o$vol ^ 2 / 2) * (t1)) / (o$vol * sqrt(t1))
  d1 = d2 + o$vol* sqrt(t1)
  d2n = ((log(o$S0/o$K) + (o$r-o$q) * t1 - o$vol^ 2 * t2 / 2) /(o$vol * sqrt(t2)))
  d1n = d2n + o$vol* sqrt(t2)
  o$PxBS =(o$S0 * exp (-o$q * t1) * stats::pnorm(d1) - o$K * exp(-o$r * t1)* stats::pnorm(d2)  
           + o$K * exp(-o$r * t1) * stats::pnorm(-d2n) - o$S0 * exp (-o$q *t1) * stats::pnorm(-d1n))
  return(o)
}


#' @title Chooser option valuation via Monte Carlo (MC) simulations
#' @description Price chooser option using Monte Carlo (MC) simulation.
#'
#' @author Xinnan Lu, Department of Statistics, Rice University, Spring 2015
#' 
#' @param o An object of class \code{OptPx}
#' @param T1 The time when the choice is made whether the option is a call or put
#' @param NPaths The number of Monte Carol simulation paths
#' @param isEu Values \code{TRUE} or \code{FALSE} indicating 
#'  if the chooser is an European or American style option
#' @param plot Values \code{TRUE} or \code{FALSE} indicating 
#'   whether to include a comparison plot of option price versus number of paths
#' @return A list of class \code{ChooserMC} consisting of original \code{OptPx} object, 
#' option pricing parameters \code{isEu}, \code{NPaths}, and \code{T1}, 
#' as well as the computed price \code{PxMC} for the chooser option.
#' 
#' 
#' @references Hull, John C., \emph{Options, Futures and Other Derivatives}, 9ed, 2014. Prentice Hall. 
#' ISBN 978-0-13-345631-8, \url{http://www-2.rotman.utoronto.ca/~hull/ofod/index.html}.
#' p.603.
#' 
#' @details 
#' A chooser option (sometimes referred to as an as you like it option) has the feature that, 
#' after a specified period of time, the holder can choose whether the option is a call or a put.
#' In this algorithm, we can price chooser options when the underlying options are both European or are both American.
#' When the underlying is an American option, the option holder can exercise before and after T1.
#' 
#' @examples
#' (o = ChooserMC())$PxMC
#' 
#' o = OptPx(Opt(Right='Call',Style="Chooser"))
#'  ChooserMC(o,isEu=TRUE,NPaths=5, plot=TRUE)
#'  
#'  o = OptPx(Opt(Right='Put',Style="Chooser"))
#'  ChooserMC(o,isEu=TRUE,NPaths=5, plot=TRUE)
#'  
#'  o = Opt(Right='C',S0=100,K=110,ttm=4,Style="Chooser")
#'  o = OptPx(o,vol=0.2,r=0.05,q=0.04)
#'  ChooserMC(o,isEu=TRUE,T1=2,NPaths=5)
#'  
#'  o = Opt(Right='P',S0=110,K=100,ttm=4,Style="Chooser")
#'  o = OptPx(o,vol=0.2,r=0.05,q=0.04)
#'  ChooserMC(o,isEu=TRUE,T1=2,NPaths=5)
#'  
#'  o = Opt(Right='C',S0=50,K=50,ttm=0.5,Style="Ch")
#'  o = OptPx(o,vol=0.25,r=0.08,q=0.1)
#'  ChooserMC(o,isEu=FALSE,T1=0.25,NPaths=5)
#'  
#' @export
#'
ChooserMC = function(o=OptPx(Opt(Style='Chooser')), isEu=TRUE, T1=1, NPaths=5, plot=FALSE){
  
  stopifnot(is.OptPx(o), o$Style$Chooser, is.numeric(T1), is.numeric(NPaths), T1<=o$ttm)
  
  ttm=o$ttm;r=o$r;q=o$q;vol=o$vol;S0=o$S0;ttm=o$ttm;K=o$K;SignCP=o$Right$SignCP;
  
  ts=(1:(252*ttm))/252 # all the trading days before maturity
  
  # simulate underlying stock prices
  Sts.NPaths = sapply(1:NPaths, function(x) sapply(ts, function(t) S0*exp((r-vol^2/2)*t+vol*stats::rnorm(1,0,1)*sqrt(t))))
  
  Max0 = function(x){1/2*abs(x)+1/2*x}
  
  PxFun = function(Sts){
    if(isEu){## European option
      # max(c,NPaths): c is the payoff of the call option at maturity (last price), discounting back to origin
      Pxs = (Max0(Sts[nrow(Sts),]-K)+exp(-q*(ttm-T1))*Max0(K*exp(-(r-q)*(ttm-T1))-S0))*exp(-(r-q)*T1)
      # find the MC simulation average option price
      Px = mean(Pxs)
    } else {## American Option
      # before the choice being made
      Pxs1 = sapply(1:(T1*252), function(x) Max0(SignCP*(Sts[x,]-K))*exp(-(r-q)*x/252))
      # after the choice being made
      # max(c,NPaths): c is all the payoffs of the call option at all the prices (Sts) before maturity
      Pxs2 = sapply((T1*252+1):nrow(Sts), function(x) Max0(Sts[x,]-K)+exp(-q*(x/252-T1))*Max0(K*exp(-(r-q)*(x/252-T1))-S0))*exp(-(r-q)*T1)
      # combine two chunks of prices
      Pxs = cbind(Pxs1,Pxs2)
      # find the early exercise where the value of the option is the maximum
      Pxs.max = apply(Pxs, 1, max) 
      Px = mean(Pxs.max)
    }
    return(Px)
  }
  
  # calculate option price paths
  PxPaths = sapply(2:NPaths, function(n) PxFun(Sts.NPaths[,1:n]))
  
  o$isEu=isEu; o$isAm=!isEu; o$NPaths=NPaths; o$T1=T1; o$PxMC=PxPaths[length(PxPaths)];
  
  # plot option prices vs. number of steps
  if (plot){
    graphics::plot(2:NPaths, PxPaths, type="n",
         main="Chooser option price convergence", xlab="Number of Paths", ylab="Price (MC)")
    graphics::lines(2:NPaths,PxPaths)
  }
  
  return(o)
  
}



#' @title Chooser option valuation via Lattice Tree (LT) Model
#' @description Calculates the price of a Chooser option using a recombining binomial tree model. 
#' Has pricing capabilities for both simple European Chooser options 
#' as well as American Chooser Options, where exercise can occur any time as a call or put options.
#'
#' @author Richard Huang, Department of Statistics, Rice University, spring 2015
#' 
#' @param o The \code{OptPx} option object to price. 
#' @param t1  The time to maturity of the call option, measured in years. 
#' @param t2  The time to maturity of the put option, measured in years. 
#' @param IncBT \code{TRUE/FALSE} Choice of including the lattice tree simulation in the output. 
#' Input \code{FALSE} yields faster computation and fewer calculated results to store in memory.
#' 
#' @details The American chooser option is interpreted as exercise of option being available 
#' at any point in time during the life of the option.
#' 
#' @return An original \code{OptPx} object with \code{PxLT} field as the price of the option and user-supplied \code{ttc}, 
#' \code{IncBT} parameters attached.
#'
#' @references Hull, J.C., \emph{Options, Futures and Other Derivatives}, 9ed, 2014. Prentice Hall. 
#' ISBN 978-0-13-345631-8, \url{http://www-2.rotman.utoronto.ca/~hull/ofod/index.html}
#' \cr Thomas S.Y. Ho et al., \emph{The Oxford Guide to Financial Modeling : Applications for Capital Markets. . .}
#'
#' @examples
#' (o = ChooserLT())$PxLT    #Default Chooser option price. (See Ho pg 234 in references)
#' 
#' o = Opt('Eu', S0=100, ttm=1, K=100)
#' o = OptPx(o, r=0.10, q=0, vol=0.1, NSteps=5)
#' (o = ChooserLT(o, t1 = .5, t2 =.5, IncBT=TRUE))$PxLT
#' 
#' #American Chooser, higher price than European equivalent
#' o = Opt('Am', S0=100, ttm=1, K=100)
#' o = OptPx(o, r=0.10, q=0, vol=0.1, NSteps=5)
#' ChooserLT(o,t1=.5, t2=.5,IncBT=FALSE)$PxLT
#' 
#' o = Opt('Eu', S0=50, ttm=1, K=50)
#' o = OptPx(o, r=0.05, q=0.02, vol=0.25, NSteps=5)
#' ChooserLT(o, t1 = .75, t2 = .75, IncBT=FALSE)$PxLT 
#' 
#' o = Opt('Eu', S0=50, ttm=1, K=50)
#' o = OptPx(o, r=0.05, q=0.5, vol=0.25, NSteps=5)
#' ChooserLT(o, t1 = .75, t2 = .75, IncBT=FALSE)$PxLT
#'  
#' @export
#' 
ChooserLT = function(
  o=OptPx(Opt('Chooser', ttm = 1)), t1 = .5, t2 = .5, IncBT=FALSE)
  { #o=OptChooserLT()
  
  stopifnot(is.OptPx(o), is.numeric(t1), is.numeric(t2), t1 == t2); # algorithm requires that a OptChooserLT object is provided
  NSteps=ceiling(t1/o$ttm * o$NSteps)  #(t1/o$ttm * o$NSteps)%%1 == 0
  p=o$p; K=o$K
  
  step.ttc = t1/o$ttm * NSteps
  S = with(o, S0*d^(0:NSteps)*u^(NSteps:0)) # vector of terminal stock prices, highest to lowest (@t=ttm)
  O.call = pmax((S - K), 0) # vector of terminal call option payouts (@t=ttm)
  O.put = pmax(-(S - K), 0) # vector of terminal put option payouts (@t=ttm)
  O.chooser = c()
  
  ReCalc.O_S.on.Prior.Time.Step = function(i) { 
    if (i >= step.ttc + 1) {
      O.call <<- o$DF_dt * (p*O.call[-i-1] + (1-p)*O.call[-1])  #prior call option prices (@time step=i-1)
      O.put <<- o$DF_dt * (p*O.put[-i-1] + (1-p)*O.put[-1])     #prior put option prices (@time step =i-1)
      S <<- o$d * S[-i-1]                   # prior stock prices (@time step=i-1)
      
      Payout.call = pmax((S - K), 0)   # call payout at time step i-1 (moving backward in time)
      Payout.put = pmax(-(S - K), 0)   # put payout at time step i-1 (moving backward in time)
      if (o$Style$American == TRUE) { 
        O.call <<- pmax(O.call, Payout.call)    #Test for American Chooser, max of payout values & expected values 
        O.put <<- pmax(O.put, Payout.put)
      }
      return(cbind(S, O.call, O.put))
    }
    else if (i == step.ttc) { # Time step where choice between put/call is made
      S <<- o$d * S[-i-1]                   # prior stock prices (@time step=i-1)
      O.chooser <<- o$DF_dt * (p*pmax(O.call[-i-1], O.put[-i-1]) + (1-p)*pmax(O.call[-1], O.put[-1]))  #prior chooser option prices (@time step=i-1)        
      Payout.call = pmax((S - K), 0)   
      Payout.put = pmax(-(S - K), 0)
      if (o$Style$American == TRUE) {
        O.chooser <<- pmax(O.chooser, Payout.call, Payout.put)    #
      }
      return(cbind(S,O.chooser))
    }
    else { #Rolling back to current time
      S <<- o$d * S[-i-1]                   
      O.chooser <<- o$DF_dt * (p*O.chooser[-i-1] + (1-p)*O.chooser[-1])  
      Payout.call = pmax((S - K), 0)  
      Payout.put = pmax(-(S - K), 0)
      if (o$Style$American == TRUE) {
        O.chooser <<- pmax(O.chooser, Payout.call, Payout.put)    #
      }
      return(cbind(S,O.chooser))
    }
  }
  
  BT = append(list(cbind(S, O.call, O.put)), sapply(NSteps:1, ReCalc.O_S.on.Prior.Time.Step)) #binomial tree
  
  #Return
  o$t1 <- t1
  o$t2 <- t2
  o$IncBT <- IncBT
  o$PxLT = BT[[length(BT)]][[2]]  # add ChoosterLT price
  if (IncBT) o$BT = BT
  return(o) 
}
