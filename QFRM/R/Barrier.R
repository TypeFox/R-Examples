#' @title Barrier option pricing via Black-Scholes (BS) model 
#' @description This function calculates the price of a Barrier option. This price is based on the assumptions that
#' the probability distribution is lognormal and that the asset price is observed continuously.
#'  @details To price the barrier option, we need to know whether the option is Up or Down | In or Out | Call or Put. Beyond that 
#'  we also need the S0, K, r, q, vol, H, and ttm arguments from the object classes defined in the package.
#'   
#' 
#'  @author Kiryl Novikau, Department of Statistics, Rice University, Spring 2015
#'  @param o The \code{OptPx} option object to price. See \code{OptBarrier()}, \code{OptPx()}, and \code{Opt()} for more information.
#'  @param dir The direction of the option to price. Either Up or Down.
#'  @param knock Whether the option goes In or Out when the barrier is reached.
#'  @param H The barrier level
#'  @return The price of the barrier option \code{o}, which is based on the BSM-adjusted algorithm (see references). 
#'  
#' @references Hull, John C., \emph{Options, Futures and Other Derivatives}, 9ed, 2014. Prentice Hall. 
#' ISBN 978-0-13-345631-8. \url{http://www-2.rotman.utoronto.ca/~hull/ofod/index.html}.
#'  pp.606-607
#' 
#'  @examples  
#'  (o = BarrierBS())$PxBS # Option with default arguments is valued at $9.71
#'  
#'  #Down-and-In-Call
#'  o = Opt(Style='Barrier', S0=50, K=50, ttm=1, Right="Call", ContrSize=10)
#'  o = OptPx(o,  r = .05, q = 0, vol = .25)
#'  o = BarrierBS(o, dir = "Down", knock = 'In', H = 40)
#'  
#'  #Down-and-Out Call
#'  o = Opt(Style='Barrier', S0=50, K=50, ttm=1, Right="Call", ContrSize=10)
#'  o = OptPx(o, r = .05, q = .02, vol = .25)
#'  o = BarrierBS(o, dir = "Down", knock = 'Out', H = 40)
#'  
#'  #Up-and-In Call
#'  o = Opt(Style='Barrier', S0=50, K=50, ttm=1, Right="Call", ContrSize=1)
#'  o = OptPx(o, r = .05, q = .02, vol = .25)
#'  o = BarrierBS(o, dir = "Up", knock = 'In', H = 60)
#'  
#'  #Up-and-Out Call
#'  o = Opt(Style='Barrier', S0 = 50, K = 50, ttm = 1, Right="Call", ContrSize=1)
#'  o = OptPx(o, r = .05, q = .02, vol = .25)
#'  o = BarrierBS(o, dir = "Up", knock = 'Out', H = 60)
#'
#'  #Down-and-In Put
#'  o = Opt(Style='Barrier', S0=50, K=50, ttm=1, Right="Put", ContrSize=1)
#'  o = OptPx(o, r = .05, q = .02, vol = .25)
#'  o = BarrierBS(o, dir = "Down", knock = 'In', H = 40)
#'  
#'  #Down-and-Out Put
#'  o = Opt(Style='Barrier', S0=50, K=50, ttm=1, Right="Put", ContrSize=1)
#'  o = OptPx(o, r = .05, q = .02, vol = .25)
#'  o = BarrierBS(o, dir = "Down", knock = 'Out', H = 40)
#'  
#'  #Up-and-In Put
#'  o = Opt(Style='Barrier', S0=50, K=50, ttm=1, Right="Put", ContrSize=1)
#'  o = OptPx(o, r = .05, q = .02, vol = .25)
#'  o = BarrierBS(o, dir = "Up", knock = 'In', H = 60)
#'  
#'  #Up-and-Out Put
#'  o = Opt(Style='Barrier', S0=50, K=50, ttm=1, Right="Put", ContrSize=1)
#'  o = OptPx(o, r = .05, q = .02, vol = .25)
#'  o = BarrierBS(o, dir = "Up", knock = 'Out', H = 60)
#'  
#'  @export
#'  
BarrierBS <- function(o = OptPx(Opt(Style = "Barrier")), dir = c("Up","Down"), knock = c('In','Out'), H = 40){
  stopifnot(o$Style$Barrier, is.OptPx(o), is.numeric(H), is.character(dir), is.character(knock))  
  o <- BS(o)
  d1 <- o$BS$d1
  d2 <- o$BS$d2
  dir=match.arg(dir)
  knock = match.arg(knock)
  isIn <- switch(knock, In = TRUE, Out = FALSE)
  isUp <- switch(dir, Up = TRUE, Down = FALSE)
  o = append(o, list(H = H, dir = dir, knock = knock, isIn = isIn, isUp = isUp, isOut = !isIn, isDown = !isUp))
  lambda <- (o$r - o$q + ((o$vol)^2 / 2)) / (o$vol)^2
  y <- (log((o$H^2)/(o$S0 * o$K)) / (o$vol * sqrt(o$ttm))) + (lambda * o$vol * sqrt(o$ttm))
  x1 <- (log(o$S0 / o$H) / (o$vol * sqrt(o$ttm))) + (lambda * o$vol * sqrt(o$ttm))
  y1 <- (log(o$H / o$S0) / (o$vol * sqrt(o$ttm))) + (lambda * o$vol * sqrt(o$ttm))        
  if(o$Right$Call){
    c = o$S0 * exp(-o$q * o$ttm) * stats::pnorm(d1) - o$K * exp(-o$r * o$ttm) * stats::pnorm(d2)
    if(o$isDown){
      if(o$H <= o$K){
        c_di <- (o$S0 * exp(-1 * o$q * o$ttm) * (o$H / o$S0)^(2*lambda) * stats::pnorm(y)) 
        - (o$K * exp(-1*o$r * o$ttm) * (o$H / o$S0)^(2 * lambda - 2) * stats::pnorm(y - (o$vol * sqrt(o$ttm))))
        if(o$isIn){
          o$PxBS <- c_di        
          return(o)
        }
        else{
          o$PxBS <- c-c_di
          return(o)
        }
      }
      if(o$H >= o$K){
        c_do <- (o$S0 * stats::pnorm(x1) * exp(-o$q * o$ttm)) - (o$K * exp(-o$r * o$ttm) * stats::pnorm(x1 - o$vol * sqrt(o$ttm))) - 
          o$S0 * exp(-o$q * o$ttm) * (o$H / o$S0)^(2*lambda) * stats::pnorm(y1) 
        + o$K*exp(-o$r * o$ttm)*(o$H/o$S0)^(2*lambda - 2) * stats::pnorm(y1 - o$vol * sqrt(o$ttm))
        if(o$isOut){
          o$PxBS <- c_do
          return(o)
        }
        else{
          o$PxBS <- c - c_do
          return(o)
        }
      }
    }
    if(o$isUp){
      if(o$H <= o$K){
        if(o$isOut){
          o$PxBS <- 0
          return(o)
        }
        else{
          o$PxBS <- c
          return(o)
        }
      }
      if(o$H >= o$K){
        c_ui <- o$S0 * stats::pnorm(x1) * exp(-o$q * o$ttm) - o$K * exp(-o$r * o$ttm) * stats::pnorm(x1 - o$vol * sqrt(o$ttm)) 
        - o$S0 * exp(-o$q * o$ttm) * (o$H / o$S0)^(2*lambda) * 
          (stats::pnorm(-y) - stats::pnorm(-y1)) 
        + o$K * exp(-o$r * o$ttm) * (o$H / o$S0)^(2*lambda - 2) * (stats::pnorm(-y + o$vol * sqrt(o$ttm)) 
                                                                   - stats::pnorm(-y1 + o$vol * sqrt(o$ttm)))
        if(o$isIn){
          o$PxBS <- c_ui
          return(o)
        }
        if(o$isOut){
          c_uo <- c - c_ui
          o$PxBS <- c_uo
          return(o)
        }
      }
      
    }
  }
  if(o$Right$Put){
    p = o$K * exp(-o$r * o$ttm) * stats::pnorm(-d2) - o$S0 * exp(-o$q * o$ttm) * stats::pnorm(-d1)
    if(o$isUp){
      if(o$H >= o$K){
        p_ui = -o$S0 * exp(-o$q * o$ttm) * (o$H / o$S0)^(2*lambda) * stats::pnorm(-y) 
        + o$K*exp(-o$r * o$ttm) * (o$H / o$S0)^(2*lambda - 2) * stats::pnorm(-y + o$vol * sqrt(o$ttm))
        if(o$isIn){
          o$PxBS <- p_ui
          return(o)
        }
        if(o$isOut){
          o$PxBS <- p - p_ui
          return(o)
        }
      }
      if(o$H <= o$K){
        p_uo = -o$S0 * stats::pnorm(-x1) * exp(-o$q * o$ttm) 
        + o$K * exp(-o$r * o$ttm) * stats::pnorm(-x1 + o$vol * sqrt(o$ttm)) +
          o$S0 * exp(-o$q * o$ttm) * (o$H / o$S0)^(2*lambda) * stats::pnorm(-y1) 
        - o$K*exp(-o$r * o$ttm) *(o$H / o$S0)^(2*lambda-2) * stats::pnorm(-y1 + o$vol * sqrt(o$ttm))
        if(o$isOut){
          o$PxBS <- p_uo
          return(o)
        }
        if(o$isIn){
          o$PxBS <- p - p_uo
          return(o)
        }
      }
    }
    if(o$isDown){
      if(o$H >= o$K){
        if(o$isIn){
          o$PxBS <- p
          return(o)
        }
        if(o$isOut){
          o$PxBS <- 0
          return(o)
        }
      }
      if(o$H <= o$K){
        p_di <- (-o$S0 * stats::pnorm(-x1) * exp(-o$q * o$ttm)) + (o$K * exp(-o$r * o$ttm) * stats::pnorm(-x1 + o$vol * sqrt(o$ttm))) +
          (o$S0 * exp(-o$q * o$ttm) * (o$H / o$S0)^(2*lambda) * (stats::pnorm(y) - stats::pnorm(y1))) -
          (o$K * exp(-o$r * o$ttm) * (o$H / o$S0)^(2*lambda - 2) * (stats::pnorm(y - o$vol * sqrt(o$ttm)) - stats::pnorm(y1 - o$vol * sqrt(o$ttm))))
        if(o$isIn){
          o$PxBS <- p_di
          return(o)
        }
        if(o$isOut){
          o$PxBS <- p - p_di
          return(o)
        }
      }
    }
  }
}



#' @title Barrier option valuation via Monte Carlo (MC) simulation.
#' @description Calculates the price of a Barrier Option using 10000 Monte Carlo simulations.
#' The helper function BarrierCal() aims to calculate expected payout for each stock prices.
#' 
#' Important Assumptions: 
#' The option follows a General Brownian Motion (GBM) 
#' \eqn{ds = mu * S * dt + sqrt(vol) * S * dW} where \eqn{dW ~ N(0,1)}. 
#' The value of \eqn{mu} (the expected percent price increase) is assumed to be \code{o$r-o$q}. 
#' 
#' @author Huang Jiayao, Risk Management and Business Intelligence at Hong Kong University of Science and Technology, 
#'          Exchange student at Rice University, Spring 2015
#' @param o The \code{OptPx} Barrier option to price.
#' @param knock Defines the Barrier option to be "\code{In}" or "\code{Out}" 
#' @param B The Barrier price level
#' @param NPaths The number of simulation paths to use in calculating the price
#' @return The option \code{o} with the price in the field \code{PxMC} based on MC simulations and the Barrier option 
#'         properties set by the users themselves 
#'  
#' @references Hull, John C., \emph{Options, Futures and Other Derivatives}, 9ed, 2014. Prentice Hall. 
#' ISBN 978-0-13-345631-8, \url{http://www-2.rotman.utoronto.ca/~hull/ofod/index.html}.
#' Also, 
#' \url{http://stackoverflow.com/questions/25946852/r-monte-carlo-simulation-price-path-converging-volatility-issue}
#'              
#' @examples
#'  (o = BarrierMC())$PxMC #Price =~ $11
#'   
#'  o = OptPx(o=Opt(Style='Barrier'),NSteps = 10)
#'  (o = BarrierMC(o))$PxMC #Price =~ $14.1
#'   
#'  (o = BarrierMC(NPaths = 5))$PxMC # Price =~ $11
#'   
#'  (o = BarrierMC(B=65))$PxMC # Price =~ $10
#'  
#'  (o = BarrierMC(knock="Out"))$PxMC #Price =~ $1
#' @export

BarrierMC = function(o=OptPx(o=Opt(Style='Barrier')),knock = c("In","Out"),B=60,NPaths=5){
  
  #First check if the inputs are valid
  stopifnot(o$Style$Barrier, is.OptPx(o), is.numeric(B),is.character(knock))
  
  K=o$K
  S0=o$S0
  
  Knock = match.arg(knock)
  isIn = switch(Knock,In=TRUE,Out=FALSE)
  
  # Decide the barrier option is up or down
  if (S0<=B) isUp=TRUE else isUp=FALSE
  
  # This helper function calculates Barrier payout given a stock price
  BarrierCal = function(x){
    
    if(isIn){
      if(isUp) {if (x>=B) payout = max(x-K,0) else payout=0  }
      else {if (x<=B) payout = max(x-K,0) else payout=0  }
    }
    else {
      if(isUp) {if (x>=B) payout = 0 else payout= max(x-K,0)}
      else {if (x<=B) Payout = 0 else payout= max(x-K,0)}
    }
    
    return(payout)
  }
  
  o$Knock = Knock
  o$B = B
  
  # Calculate Barrier option price using Monte Carlo Simulations
  
  o$PxMC = o$Right$SignCP*mean(sapply((1:NPaths),function(trial_sum){
    div = with(o, exp(((r-q) - 0.5 * vol^2) * dt + vol * sqrt(dt) * rnorm(NSteps)))
    inprices = cumprod(c(o$S0, div))
    
    payout = max(sapply(inprices,BarrierCal))
    return(exp(-o$r * o$ttm) * payout)
  }))
  
  return(o)
}



#' @title Barrrier option valuation via lattice tree (LT)
#' @description Use Binomial Tree to price barrier options with relatively large NSteps (NSteps > 100) steps.
#'              The price may be not as percise as BSM function cause the convergence speed for Binomial Tree is kind of slow.
#' @author Tong Liu, Department of Statistics, Rice University, Spring 2015
#' 
#' @param o      An object of class \code{OptPx}
#' @param dir    A direction for the barrier, either \code{'Up'} or \code{'Down'} Default=\code{'Up'}
#' @param knock  The option is either a knock-in option or knock-out option. Default=\code{'In'}
#' @param H      The barrier level. \code{H} should less than \code{S0} if \code{'Up'}, 
#'   \code{H} should greater than \code{S0} if \code{'Down'} Default=60.
#' @return A list of class \code{BarrierLT} consisting of the input object \code{OptPx} 
#' and the appended new parameters and option price.
#'  
#'  @references Hull, John C., \emph{Options, Futures and Other Derivatives}, 9ed, 2014. Prentice Hall. 
#'  ISBN 978-0-13-345631-8, \url{http://www-2.rotman.utoronto.ca/~hull/ofod/index.html}
#' \cr p.467-468. Trinomial Trees, p.604-606: Barrier Options.
#' @examples 
#' # default Up and Knock-in Call Option with H=60, approximately 7.09
#' (o = BarrierLT())$PxLT   
#' 
#' #Visualization of price changes as Nsteps change.
#' o = Opt(Style="Barrier")
#' visual=sapply(10:200,function(n) BarrierLT(OptPx(o,NSteps=n))$PxLT)
#'  
#' c=(10:200)
#' plot(visual~c,type="l",xlab="NSteps",ylab="Price",main="Price converence with NSteps")
#' 
#' # Down and Knock-out Call Option with H=40
#' o = OptPx(o=Opt(Style="Barrier"))
#' BarrierLT(o,dir="Down",knock="Out",H=40)  
#' 
#' # Down and Knock-in Call Option with H=40
#' o = OptPx(o=Opt(Style="Barrier"))
#' BarrierLT(o,dir="Down",knock="In",H=40) 
#' 
#' # Up and Knock-out Call Option with H=60
#' o = OptPx(o=Opt(Style="Barrier"))
#' BarrierLT(o,dir='Up',knock="Out") 
#' 
#' # Down and Knock-out Put Option with H=40
#' o = OptPx(o=Opt(Style="Barrier",Right="Put"))
#' BarrierLT(o,dir="Down",knock="Out",H=40)
#' 
#' # Down and Knock-in Put Option with H=40  
#' o = OptPx(o=Opt(Style="Barrier",Right="Put"))
#' BarrierLT(o,dir="Down",knock="In",H=40) 
#' 
#' # Up and Knock-out Put Option with H=60
#' o = OptPx(o=Opt(Style="Barrier",Right="Put"))
#' BarrierLT(o,dir='Up',knock="Out")
#' 
#' # Up and Knock-in Put Option with H=60
#' BarrierLT(OptPx(o=Opt(Style="Barrier",Right="Put"))) 
#' @export
#' 
BarrierLT=function(o=OptPx(Opt(Style="Barrier"),vol=0.25,r=0.05,q=0.02,NSteps=5),
                   dir = c('Up', 'Down'), knock = c('In', 'Out'), H = 60){
  #Condition 
  stopifnot(is.OptPx(o), o$Style$Name=='Barrier',is.character(dir), is.character(knock), is.numeric(H))
  #---Direction and Knock
  dir=match.arg(dir)
  Knock=match.arg(knock)
  is.Up=switch(dir,Up=TRUE,Down=FALSE)
  is.In=switch(Knock,In=TRUE,Out=FALSE)
  is.Out=!is.In
  is.Down=!is.Up
  
  o.class = class(o)    # remember the class name
  
  # Regular Options Price
  st <- pmax (o$S0*o$d^(o$NSteps:0)*o$u^(0:o$NSteps)-o$K,0)
  pt <- stats::dbinom(0:o$NSteps,o$NSteps,o$p)
  Px.Vanilla <- o$DF_ttm*sum(st*pt)
  
  # First Calculate Knock-out, then use parity to calculate Knock-in.
  S = with(o, o$S0*o$d^(0:o$NSteps)*o$u^(o$NSteps:0)) # vector of terminal stock prices, lowest to highest (@t=ttm)
  if (is.Down)   O=(S>H)*pmax(o$Right$SignCP*(S-o$K),0) # vector of terminal option value
  else           O=(S<H)*pmax(o$Right$SignCP*(S-o$K),0)
  
  ReCalc.O_S.on.Prior.Time.Step = function(i) { #sapply(1:(i-1), function(j) 
    S <<- S[1:i]/o$u  #vector of prior stock prices (@time step=i-1)
    if (is.Up){
      stopifnot(o$S0<H)         #if Up, S0 should smaller than H, otherwise meaningless
      O <<- (S<H)*(O[1:i]*o$p+O[2:(i+1)]*(1-o$p))*o$DF_dt
    }
    else if (is.Down){
      stopifnot(o$S0>H)         #if Down, S0 should greater than H, otherwise meaningless
      O <<- (S>H)*(O[1:i]*o$p+O[2:(i+1)]*(1-o$p))*o$DF_dt
    }
    return(cbind(S, O))
  }
  
  BT=append(list(cbind(S, O)), sapply(o$NSteps:1, ReCalc.O_S.on.Prior.Time.Step,simplify=F))
  Px.Out = BT[[length(BT)]][[2]]  # add BOPM price
  
  o$Direction=dir
  o$Knock=Knock
  o$isUp=is.Up
  o$isIn=is.In
  
  if (is.Out)       o$PxLT=Px.Out
  else              o$PxLT=Px.Vanilla-Px.Out
  
  class(o)=o.class
  return(o) 
}


