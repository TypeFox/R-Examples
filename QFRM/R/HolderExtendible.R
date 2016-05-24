#'@title        Holder Extendible option valuation via Black-Scholes (BS) model
#'@description  Computes the price of exotic option (via BS model)
#'  which gives the holder the right to extend the option's maturity at an additional premium.
#'
#'@author       Le You, Department of Statistics, Rice University, Spring 2015
#'
#'@param        o An object of class \code{OptPx}
#'@param        t1    The time to maturity of the call option, measured in years.
#'@param        t2    The time to maturity of the put option, measured in years.
#'@param        k     The exercise price of the option at t2, a numeric value.
#'@param        A     The corresponding asset price has exceeded the exercise price X.
#'@return       The original \code{OptPx} object 
#'  and the option pricing parameters \code{t1}, \code{t2},\code{k},\code{A}, and computed price \code{PxBS}.
#'
#'@references   Hull, J.C., \emph{Options, Futures and Other Derivatives}, 9ed, 2014. Prentice Hall. ISBN 978-0-13-345631-8, 
#'              \url{http://www-2.rotman.utoronto.ca/~hull/ofod/index.html}
#'              \cr Haug, Espen G.,\emph{Option Pricing Formulas}, 2ed.
#'              
#'@examples
#' (o = HolderExtendibleBS())$PxBS
#' 
#' o = Opt(Style='HolderExtendible',Right='Call', S0=100, ttm=0.5, K=100)
#' o = OptPx(o,r=0.08,q=0,vol=0.25)
#' (o = HolderExtendibleBS(o,k=105,t1=0.5,t2=0.75,A=1))$PxBS
#' 
#' o = Opt("HolderExtendible","Put", S0=100, ttm=0.5, K=100)
#' o = OptPx(o,r=0.08,q=0,vol=0.25)
#' (o = HolderExtendibleBS(o,k=90,t1=0.5,t2=0.75,A=1))$PxBS
#'
#' @export
#'
HolderExtendibleBS = function(o=OptPx(Opt(Style='HolderExtendible')),k=105,t1=0.5,t2=0.75,A=1)
  { 
    isCall=as.numeric(o$Right$Call)
    isPut=as.numeric(o$Right$Put)  
    I1 = stats::uniroot(function(I1) BS(OptPx(Opt(Style='HolderExtendible',Right=o$Right$Name,S0=I1,K=k,ttm=abs(t2-t1)),r=o$r,q=o$q,vol=o$vol))$PxBS- A-isPut*o$K+isPut*I1, 
                 interval=c(0,200), tol=0.000001)$root
    I2 = stats::uniroot(function(I2) BS(OptPx(Opt(Style='HolderExtendible',Right=o$Right$Name,S0=I2,K=k,ttm=abs(t2-t1)),r=o$r,q=o$q,vol=o$vol))$PxBS- A+isCall*o$K-isCall*I2,
                 interval=c(0,200), tol=0.000001)$root
    y1 = (log(o$S0/I2) + ((o$r-o$q) + o$vol ^ 2 / 2) * (t1)) / (o$vol * sqrt(t1))
    y2 = (log(o$S0/I1) + ((o$r-o$q) + o$vol ^ 2 / 2) * (t1)) / (o$vol * sqrt(t1))
    
    z1 = (log(o$S0/k) + ((o$r-o$q) + o$vol ^ 2 / 2) * (t2)) / (o$vol * sqrt(t2))
    z2 = (log(o$S0/o$K) + ((o$r-o$q) + o$vol ^ 2 / 2) * (t1)) / (o$vol * sqrt(t1))
    rho = sqrt(t1/t2)
    
    o$PxBS= BS(OptPx(Opt(Style='HolderExtendible',Right=o$Right$Name,S0=o$S0,K=o$K,ttm=t1),r=o$r,q=o$q,vol=o$vol))$PxBS+o$Right$SignCP*o$S0*exp(-o$q*t2)*(pbnorm(y2,o$Right$SignCP*z1,rho)-pbnorm(y1,o$Right$SignCP*z1,rho))
    -o$Right$SignCP*k*exp(-o$r*t2)*(pbnorm(y2-o$vol*sqrt(t1),o$Right$SignCP*z1
                                           -o$Right$SignCP*o$vol*sqrt(t2),rho)
                                    -pbnorm(y1-o$vol*sqrt(t1),o$Right$SignCP*z1-o$Right$SignCP*o$vol*sqrt(t2),rho))
    -o$Right$SignCP*o$S0*exp(-o$q*t1)*(stats::pnorm(z2)-isCall*stats::pnorm(y1)-isPut*stats::pnorm(y2))
    +o$Right$SignCP*o$K*exp(-o$r*t1)*(isPut*stats::pnorm(y2-o$vol*sqrt(t1))
                                      +o$Right$SignCP*stats::pnorm(z2-o$vol*sqrt(t1))-isCall*stats::pnorm(y1-o$vol*sqrt(t1)))
    -A*exp(-o$r*t1)*(stats::pnorm(y2-o$vol*sqrt(t1))-stats::pnorm(y1-o$vol*sqrt(t1))) 
    return(o)
  }
