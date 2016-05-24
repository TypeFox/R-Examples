#' @title Shout option valuation via Monte Carlo (MC) simulations.
#' @description Calculates the price of a shout option using Monte Carlo simulations to 
#' determine expected payout. Assumes that the option follows a General
#' Brownian Motion  (GBM) process, \eqn{ds = mu * S * dt + sqrt(vol) * S * dW} where \eqn{dW ~ N(0,1)}.
#' Note that the value of \eqn{mu} (the expected price increase) is assumped to be 
#' \code{o$r}, the risk free rate of return. 
#' 
#'  @author Jake Kornblau, Department of Statistics, Rice University, 2015
#'  @param o The \code{OptPx} Shout option to price. 
#'  @param NPaths The number of simulation paths to use in calculating the price; must be >= 10
#'  
#'  @return The option object \code{o} with the price in the field \code{PxMC} based on the MC simulations.
#' 
#'  @references Hull, J.C., \emph{Options, Futures and Other Derivatives}, 9ed, 2014. Prentice Hall. 
#'  ISBN 978-0-13-345631-8, \url{http://www-2.rotman.utoronto.ca/~hull/ofod/index.html}.
#'  \cr
#'  Also: \url{http://www.math.umn.edu/~spirn/5076/Lecture16.pdf}
#'  
#'  @examples
#'   (o = ShoutMC())$PxMC # Approximately valued at $11
#'   
#'   o = Opt(Style='Shout')
#'   (o = ShoutMC(OptPx(o, NSteps = 5)))$PxMC # Approximately valued at $18.6
#'   
#'   o = Opt(Style='Shout',S0=110,K=100,ttm=.5)
#'   o = OptPx(o, r=.05, vol=.2, q=0, NSteps = 10)
#'   (o = ShoutMC(o, NPaths = 10))$PxMC
#'   @export
#'   
ShoutMC = function(o = OptPx(o = Opt(Style='Shout')), NPaths = 10) {
  
  #Stops if initial conditions are not true
  stopifnot(is.OptPx(o), o$Style$Shout, is.numeric(NPaths), NPaths>9, o$Style$Shout);
  
  sim_prices = matrix(nrow = NPaths, ncol = o$NSteps)
  V_date = matrix(data = 0, nrow = NPaths, ncol = o$NSteps)
  
  # Function for simulating the price of a stock
  ShoutSim = function(trial_num) {
    # ds = mu * S * dt + sigma * S * dW = S * (mu * dt + sigma * dW)
    ds_div_S = with(o, exp((r - 0.5 * vol^2) * dt + vol * sqrt(dt) * rnorm(NSteps)))
    
    # ds is the product of a RV and the previous price. cumprod with S0 at
    # the beginning will accomplish this.
    prices = cumprod(c(o$S0, ds_div_S))
    prices = prices[2:(length(prices))]
    
    sim_prices[trial_num, 1:(length(prices))] <<- prices
    V_date[trial_num, o$NSteps] <<- max(prices[length(prices)] - o$K, 0)
  }
  
  # Function for determining whether the option should be exercised (shouted) at any given date
  ValueCalculation = function (date) {
    shoutPossible = sim_prices[1:NPaths, date] > o$K
    
    size = sum(shoutPossible)
    
    to_optimize = matrix(data = 1, nrow = size, ncol = 3)
    to_optimize[1:size, 2] = sim_prices[shoutPossible, date]
    to_optimize[1:size, 3] = sim_prices[shoutPossible, date] ^ 2
    
    coefficients = stats::lsfit(to_optimize, V_date[shoutPossible, date + 1], intercept = FALSE)$coef
    
    values = to_optimize %*% coefficients
    V_date[shoutPossible, date] <<- (values > V_date[shoutPossible, date]) * (sim_prices[shoutPossible, date] - o$K)
  }
  
  # Calculates the price based on shout information and prices for each trial
  CalcPrice = function (trial_num) {
    shout_dates = V_date[trial_num, 1:(o$NSteps)] > 0
    price = max(o$Right$SignCP * (sim_prices[trial_num, o$NSteps] - o$K), 0)
    if (length(shout_dates) > 0) {
      price = max(price, o$Right$SignCP *(sim_prices[trial_num, min(shout_dates)]- o$K))
    }
    
    return(exp(-o$r * o$ttm) * price)
  }
  
  o$NPaths = NPaths
  sapply((1:NPaths), ShoutSim)
  if (o$NSteps > 1) {
    sapply(((o$NSteps - 1): 1), ValueCalculation)
  }
  
  o$PxMC = mean(sapply((1:NPaths), CalcPrice))
  
  return(o)
}


#' @title       Shout option valuation via lattice tree (LT) 
#' @author      Le You, Department of Statistics, Rice University, Spring 2015
#' @description A shout option is a European option where the holder can shout to the writer at one time during its life.
#'              At the end of the life of the option, the option holder receives either the usual payoff from a European option
#'              or the instrinsic value at the time of the shout, which ever is greater. 
#'              \eqn{max(0,S_T-S_tau)+(S_tau-K)}
#' 
#' @param       o An object of class \code{OptPx}
#' @return       A list of class \code{ShoutLT} consisting of the original \code{OptPx} object, 
#'              binomial tree step \code{BT} and the computed price \code{PxBS}.
#'              
#' @references Hull, J.C., \emph{Options, Futures and Other Derivatives}, 9ed, 2014. Prentice Hall. 
#' ISBN 978-0-13-345631-8, \url{http://www-2.rotman.utoronto.ca/~hull/ofod}. \url{http://amzn.com/0133456315}
#' 
#' 
#' @examples 
#' (o = ShoutLTVectorized( OptPx(Opt(Style='Shout'))))$PxLT
#' 
#' o = Opt(Style='Shout')
#' (o = ShoutLTVectorized( OptPx(o, r=.1, q=.02, vol=.45, NSteps=10)))$PxLT
#' 
#' @export
#'
ShoutLTVectorized = function(o = OptPx(o = Opt(Style='Shout'))){
  stopifnot(is.OptPx(o), o$Style$Shout)
  
  st = pmax(o$S0*o$d^(o$NSteps:0)*o$u^(0:o$NSteps)-o$K,0)
  pt = stats::dbinom(0:o$NSteps,o$NSteps,o$p)
  o$PxLT = o$DF_ttm*sum(st*pt)
  return (o)
}



#' @title       Shout option valuation via lattice tree (LT) 
#' @description A shout option is a European option where the holder can shout to the writer at one time during its life.
#'              At the end of the life of the option, the option holder receives either the usual payoff from a European option
#'              or the instrinsic value at the time of the shout, which ever is greater.
#'              \eqn{max(0,S_T-S_tau)+(S_tau-K)}
#' @author       Le You, Department of Statistics, Rice University, Spring 2015
#'
#' @param       o An object of class \code{OptPx}
#' @param       IncBT TRUE/FALSE indicating whether to include binomial tree (list object) with output
#' @return      A list of class \code{ShoutLT} consisting of the original \code{OptPx} object,
#'              binomial tree step\code{BT} and the computed price \code{PxBS}.
#'              
#' @references Hull, J.C., \emph{Options, Futures and Other Derivatives}, 9ed, 2014. Prentice Hall. 
#' ISBN 978-0-13-345631-8, \url{http://www-2.rotman.utoronto.ca/~hull/ofod}. \url{http://amzn.com/0133456315}
#' 
#' @examples
#' (o = ShoutLT( OptPx(Opt(Style='Shout'))))$PxLT
#'
#' o = Opt(Style='Shout', Right='Call', S0=60, ttm=.25, K=60)
#' ShoutLT( OptPx(o,r=.1, q=.02, vol=.45, NSteps=10))
#'
#' o = Opt(Style='Shout', Right='Call', S0=60, ttm=.25, K=60)
#' @export
#'
ShoutLT = function(o=OptPx(Opt(Style='Shout')), IncBT=TRUE){
  stopifnot(is.OptPx(o)); # algorithm requires that a OptBT object is provided
  NSteps=o$NSteps; p=o$p; K=o$K
  if (o$Style$Shout) {
    return(ShoutLTVectorized(o)) #use vectorized calculation from above, if tree is not needed
  } else { 
    S = with(o, S0*d^(NSteps:0)*u^(0:NSteps)) # vector of terminal stock prices, lowest to highest (@t=ttm)
    O = pmax(o$SignCP * (S - K), 0) # vector of terminal option payouts (@t=ttm)
    
    ReCalc.O_S.on.Prior.Time.Step = function(i) { #sapply(1:(i-1), function(j) 
      O <<- (p*O[-i-1] + (1-p)*O[-1]) #vector of prior option prices (@time step=i-1)
      S <<- o$d * S[-i-1]  #vector of prior stock prices (@time step=i-1)
      
      Payout = pmax(o$SignCP *(S-K),0)  # payout at time step i-1 (moving backward in time)
      O <<- pmax(O, Payout)    
      return(cbind(S, O))
    }
    
    BT=append(list(cbind(S, O)), sapply(NSteps:1, ReCalc.O_S.on.Prior.Time.Step,simplify=F))
    o$PxLT = BT[[length(BT)]][[2]]  # add Shout price
    if (IncBT) o$BT = BT
    return(o)
  } 
}



#' @title Shout option valuation via finite differences (FD) method 
#' @description Shout option valuation via finite differences (FD) method
#'
#' @author Xinnan Lu, Department of Statistics, Rice University, 2015
#'
#' @param o An object of class \code{OptPx}
#' @param N The number of equally spaced intervals. Default is 100.
#' @param M The number of equally spaced stock price. Default is 20.
#' @param Smin similar to Smax
#' @param Smax A stock price sufficiently high that, when it is reached, the put option has virtually no value. 
#' The level of Smax should be chosen in such a way that one of these equally spaced stock prices is the current stock price.  
#' 
#' @references Hull, J.C., \emph{Options, Futures and Other Derivatives}, 9ed, 2014. Prentice Hall. 
#' ISBN 978-0-13-345631-8, \url{http://www-2.rotman.utoronto.ca/~hull/ofod/index.html} pp.609.
#' 
#' @details 
#' A shout option is a European option where the holder can 'shout' to the writer at one time during its life. 
#' At the end of the life of the option, the option holder receives either the usual payoff from a European option 
#' or the intrinsic value at the time of the shout, whichever is greater.
#' An explicit finite difference method (Page 482 in Hull's book) is used here to price the shout put option. 
#' Similar to pricing American options, the value of the option is consolidated at each node of the grid 
#' to see if shouting would be optimal. 
#' The corresponding shout call option is priced using the Put-Call-Parity in the finite difference method .
#' 
#' @return A list of class \code{OptPx}, including
#' option pricing parameters \code{N}, \code{M}, \code{Smin}, and \code{Smax},
#' as well as the computed option price \code{PxFD}.
#' 
#' @examples
#'  (o = ShoutFD(OptPx(Opt(Right="C", Style="Shout"))))$PxFD
#'  
#'  o = OptPx(Opt(Right="C", Style="Shout"))
#'  (o = ShoutFD(o, N=10))$PxFD # very differnt result for N=10
#'  
#'  (o = ShoutFD(OptPx(Opt(Right="P", Style="Shout"))))$PxFD
#'  
#'  o = Opt(Right='P', S0=100, K=110, ttm=0.5, Style='Shout')
#'  o = OptPx(o, vol=0.2, r=0.05, q=0.04)
#'  (o = ShoutFD(o,N=100,Smax=200))$PxFD
#'  
#'  o = Opt(Right="C", S0=110, K=100, ttm=0.5, Style="Shout")
#'  o = OptPx(o, vol=0.2, r=0.05, q=0.04)
#'  (o = ShoutFD(o,N=100,Smax=200))$PxFD
#'  
#' @export
#'
ShoutFD = function(o=OptPx(Opt(Style='Shout')), N=100, M=20, Smin=0, Smax=100){
  stopifnot (is.OptPx(o), o$Style$Shout, is.numeric(N+M+Smin+Smax))
  ttm=o$ttm; r=o$r; q=o$q; S0=o$S0; K=o$K; vol=o$vol; SignCP=o$Right$SignCP
  
  Max0 = function(x){1/2*abs(x)+1/2*x}
  
  ## discretize time and price
  dt = ttm/N
  dS = (Smax-Smin)/M
  
  ## check if S0 is j*dS and Smax sufficiently large: j = S0/dS in 1:M
  stopifnot ((S0/dS)%in%0:M, Max0(if(SignCP==1) Smin-K else K-Smax) <= 0)
  
  ## empty grid
  grid = matrix(NA, nrow=M+1, ncol=N+1, dimnames = list(S=(M:0)*dS, t=(N:0)*dt))
  ## 2 known edges of 0's and K's
  grid[rownames(grid)==Smin,] = (if(SignCP==1) 0 else K)
  grid[rownames(grid)==Smax,] = (if(SignCP==1) K else 0)
  
  ## function calc recursively filling in the grid
  calc = function(mat, n){
    # a shout vector to check if shout is optimal at each node
    shout = SignCP*(as.numeric(rownames(grid))-K)[-c(1,M+1)]
    # j vector (S) to calculate a,b,c
    j = (M-1):1 # j is from 0 to M
    
    if (n == N+1){# stop when n reaches N, the second last column
      return(sapply((M-1):1, function(x) Max0(SignCP*(x*dS-K))))
    } else {# recursion
      mat[2:M, (n+1):(N+1)] = Recall(mat, n+1)
      
      # price the option when not shout
      f = sapply(1:(M-1), function(x) {
        # explicit method 
        a = (1/(1+r*dt))*(-1/2*(r-q)*j[x]*dt+1/2*vol^2*j[x]^2*dt);
        b = (1/(1+r*dt))*(1-vol^2*j[x]^2*dt);
        c = (1/(1+r*dt))*(1/2*(r-q)*j[x]*dt+1/2*vol^2*j[x]^2*dt);
        v = t(c(a,b,c)); # vectorize
        v%*%mat[x:(x+2),n+1]})
      
      # get the maximum of the values with and without shout
      mat[2:M, n]=pmax(f, shout*exp(-r*dt*(N+1-n)))
    }  
    return(mat[2:M, n:(N+1)])
  }  
  grid[2:M, ]=calc(grid,1) # initialize recursion at n=1
  
  o$N=N; o$M=M; o$Smin=Smin; o$Smax=Smax;
  o$PxFD = grid[rownames(grid)==S0,1]  
  return(o)
}
