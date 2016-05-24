#' @title Variance Swap valuation via Black-Scholes (BS) model
#' @author Max Lee, Department of Statistics, Rice University, Spring 2015
#' 
#' @description Variance Swap valuation via Black-Scholes (BS) model
#' @param o An object of class \code{OptPx}
#' @param K A vector of non-negative strike prices
#' @param Vol a vector of non-negative, less than zero implied volatilities for the associated strikes 
#' @param notional A numeric positive amount to be invested
#' @param varrate A numeric positive varaince rate to be swapped 
#' @return An object of class \code{OptPx} with value included
#' @references Hull, J.C., \emph{Options, Futures and Other Derivatives}, 9ed, 2014. Prentice Hall. 
#' ISBN 978-0-13-345631-8, \url{http://www-2.rotman.utoronto.ca/~hull/ofod}.
#' @examples
#' (o = VarianceSwapBS())$PxBS
#' 
#' o = Opt(Style="VarianceSwap",Right="Other",ttm=.25,S0=1020)
#' o = OptPx(o,r=.04,q=.01)
#' Vol = Vol=c(.29,.28,.27,.26,.25,.24,.23,.22,.21)
#' (o = VarianceSwapBS(o,K=seq(800,1200,50),Vol=Vol,notional=10^8,varrate=.045))$PxBS
#' 
#' o = Opt(Style="VarianceSwap",Right="Other",ttm=.25,S0=1020)
#' o = OptPx(o,r=.04,q=.01)
#' Vol=c(.2,.205,.21,.215,.22,.225,.23,.235,.24)
#' (o =VarianceSwapBS(o,K=seq(800,1200,50),Vol=Vol,notional=10^8,varrate=.045))$PxBS
#' 
#' o = Opt(Style="VarianceSwap",Right="Other",ttm=.1,S0=100)
#' o = OptPx(o,r=.03,q=.02)
#' Vol=c(.2,.19,.18,.17,.16,.15,.14,.13,.12)
#' (o =VarianceSwapBS(o,K=seq(80,120,5),Vol=Vol,notional=10^4,varrate=.03))$PxBS
#' @export
#' 
VarianceSwapBS <- function(o=OptPx(Opt(Style="VarianceSwap",Right="Other",ttm=.25,S0=1020),r=.04,q=.01),K=seq(800,1200,50),Vol =seq(.2,.24,.005),notional=10^8,varrate=.045){
  #error check
  stopifnot(o$Style$VarianceSwap, is.OptPx(o),is.vector(K),is.vector(Vol),is.numeric(notional),is.numeric(varrate)) 
  
  #sort strikes and volatilities
  SANDV <- cbind(K,Vol)
  SANDV <- SANDV[order(SANDV[,1]),]
  K <- SANDV[,1]
  Vol <- SANDV[,2]
  
  #Forward Price
  F0 <- o$S0*exp((o$r-o$q)*o$ttm) 
  
  #calculate change in strikes
  deltaK <- c()  
  
  deltaK<-sapply(1:length(K),function(x){    
    if (x==1){                                 
      deltaK[x] <- K[x+1]-K[x]   
    }else if (x==length(K)){                
      deltaK[x] <- K[x]-K[x-1]             
    } else{                                    
      deltaK[x] <- .5*(K[x+1]-K[x-1])        
    } 
  })
  
  #calculate strike star
  Sstar <- c()                                  
  Sstar<-sapply(1:length(K),function(x){
    if(K[x] < F0 & K[x+1] > F0){
      Sstar <- K[x]
    }
  })
  
  Sstar <- as.vector(unlist(Sstar))  
  #calculate option values of individual strike
  optionvalues <- c()  
  
  optionvalues <- sapply(1:length(K), function(x){
    d1 <- 1/(Vol[x]*sqrt(o$ttm))*(log(o$S0/K[x])+(o$r+Vol[x]^2/2)*o$ttm)
    d2 <- d1 - Vol[x]*sqrt(o$ttm)
    
    if (K[x] < Sstar){  
      optionvalues[x] <- exp(-o$r*o$ttm)*(K[x]*stats::pnorm(-d2)-F0*stats::pnorm(-d1))
    }else if (K[x] > Sstar){
      optionvalues[x] <- exp(-o$r*o$ttm)*(F0*stats::pnorm(d1)-K[x]*stats::pnorm(d2))
    }else if (K[x] == Sstar){
      optionvalues[x] <- .5*((exp(-o$r*o$ttm)*(K[x]*stats::pnorm(-d2)-F0*stats::pnorm(-d1)))
                             +(exp(-o$r*o$ttm)*(F0*stats::pnorm(d1)-K[x]*stats::pnorm(d2))))
    }
  })
  
  integralsum <- sum(deltaK/K^2*optionvalues*exp(o$r*o$ttm)) 
  
  EV <- 2/o$ttm*log(F0/Sstar)-2/o$ttm*(F0/Sstar-1)+2/o$ttm*integralsum 
  
  VarSwap <- notional*(EV-varrate)*exp(-o$r*o$ttm) 
  o$PxBS = VarSwap
  
  return(o)
}



#' @title VarianceSwap option valuation via Monte Carlo (MC) simulation.
#' @description Calculates the price of a VarianceSwap Option using 500 Monte Carlo simulations.
#' \cr Important Assumptions: 
#' The option o followes a General Brownian Motion 
#' \eqn{ds = mu * S * dt + sqrt(vol) * S * dW} where \eqn{dW ~ N(0,1)}. 
#' The value of \eqn{mu} (the expected price increase) is assumed to be \code{o$r-o$q}. 
#' 
#' @author Huang Jiayao, Risk Management and Business Intelligence at Hong Kong University of Science and Technology, 
#'          Exchange student at Rice University, Spring 2015
#' @param o The \code{OptPx} Variance Swap option to price.
#' @param var The variance strike level
#' @param NPaths The number of simulation paths to use in calculating the price,
#' @return The option \code{o} with the price in the field \code{PxMC} based on MC simulations and the Variance Swap option 
#'         properties set by the users themselves 
#'  
#' @references Hull, J.C., \emph{Options, Futures and Other Derivatives}, 9ed, 2014. Prentice Hall. 
#' ISBN 978-0-13-345631-8, \url{http://www-2.rotman.utoronto.ca/~hull/ofod}.
#' \cr \url{http://stackoverflow.com/questions/25946852/r-monte-carlo-simulation-price-path-converging-volatility-issue}
#'              
#' @examples
#'  (o = VarianceSwapMC())$PxMC #Price = ~0.0245
#'   
#'  (o = VarianceSwapMC(NPaths = 5))$PxMC # Price = ~0.0245
#' 
#'  (o = VarianceSwapMC(var=0.4))$PxMC # Price = ~-0.1565
#' @export
#'
VarianceSwapMC=function(o=OptPx(o=Opt(Style='VarianceSwap')), var=0.2, NPaths=5){
  
  ndays = floor(252*o$ttm) #number of trading days to be simulated
  dt=1/252 
  
  o$var=var; o$NPaths = NPaths
  
  o$PxMC = o$Right$SignCP*mean(sapply((1:NPaths),function(trial_sum){
    div = with(o, exp(((r-q) - 0.5 * vol^2) * dt + vol * sqrt(dt) * rnorm(ndays)))
    inprices = cumprod(c(o$S0, div))
    
    logprices = sapply(inprices,log)
    logreturn = diff(logprices)
    logreturn2= sapply(logreturn,function(x) x^2)
    rvol = sum(logreturn2)*252*0.01/(ndays-2) #a vector of realized variance
    
    return(exp(-o$r * o$ttm) * (rvol-var))
  }))
  
  return(o)
}




