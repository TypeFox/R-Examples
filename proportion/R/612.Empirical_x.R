#' The empirical Bayesian approach for Beta-Binomial model given x
#' @param x - Number of successes
#' @param n - Number of trials 
#' @param alp - Alpha value (significance level required)
#' @param sL - Lower support for MLE optimization
#' @param sU - Upper support for MLE optimization
#' @details  Highest Probability Density (HPD) and two tailed intervals are provided for the 
#' required x (any one value from \eqn{ 0, 1, 2 ..n}) based on empirical Bayesian approach for 
#' Beta-Binomial model. Lower and Upper support values are needed to obtain the MLE of 
#' marginal likelihood for prior parameters.
#' @return A dataframe with 
#'  \item{x }{- Number of successes (positive samples)}
#'  \item{pomean }{ - Posterior mean}
#'  \item{LEBAQ }{ - Lower limits of Quantile based intervals}
#'  \item{UEBAQ }{ - Upper limits of Quantile based intervals}
#'  \item{LEBAH }{ - Lower limits of HPD intervals}
#'  \item{UEBAH }{ - Upper limits of HPD intervals}
#' @family Miscellaneous  functions for Bayesian method
#' @examples 
#' sL=runif(1,0,2)				#Lower and upper of Support for MLE optimization
#' sU=runif(1,sL,10)
#' x=0; n= 5; alp=0.05
#' empericalBAx(x,n,alp,sL,sU) 
#' @references 
#' [1] 1998 Lehmann EL and Casella G
#' Theory of Point Estimation, 2nd ed Springer, New York
#' @export
empericalBAx<-function(x,n,alp,sL,sU)
{
  if (missing(x)) stop("'x' is missing")
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(sL)) stop("'Lower support' for MLE is missing")
  if (missing(sU)) stop("'Upper support' for MLE is missing")
  if ((class(x) != "integer") & (class(x) != "numeric") || length(x) >1|| x>n || x<0 ) stop("'x' has to be between 0 and n")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<0 ) stop("'n' has to be greater than or equal to 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(sL) != "integer") & (class(sL) != "numeric") || length(sL) >1 || sL<=0 ) stop("'sL' has to be a number greater than 0")
  if ((class(sU) != "integer") & (class(sU) != "numeric") || length(sU) >1 || sU<=sL ) stop("'sU' has to be a number greater than sL")
  
    likelhd = function(y) 
    {
      a<-y[1]
      b<-y[2]
      -(choose(n,x)/beta(a,b))*beta(x+a,n-x+b)
    }
    mle=optim(c(0.1,0.1), likelhd, NULL,method = "L-BFGS-B",lower = sL, upper = sU, hessian = TRUE)$par
    
    av=mle[1]
    bv=mle[2]
    
  #  library(TeachingDemos)				#To get HPDs
    
    #####Bayesian Posterior Estimates
    #Posterior mean
    pomean=(x+av)/(n+av+bv)
    #Quantile Based Intervals
    LEBAQ=qbeta(alp/2,x+av,n-x+bv)
    UEBAQ=qbeta(1-(alp/2),x+av,n-x+bv)
    
    LEBAH=TeachingDemos::hpd(qbeta,shape1=x+av,shape2=n-x+bv,conf=1-alp)[1]
    UEBAH=TeachingDemos::hpd(qbeta,shape1=x+av,shape2=n-x+bv,conf=1-alp)[2]
    
    return(round(data.frame(x,pomean,LEBAQ,UEBAQ,LEBAH,UEBAH),4))
  }
  