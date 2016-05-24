#' The empirical Bayesian approach for Beta-Binomial model
#' @param n - Number of trials 
#' @param alp - Alpha value (significance level required)
#' @param sL - Lower support for MLE optimization
#' @param sU - Upper support for MLE optimization
#' @details  Highest Probability Density (HPD) and two tailed intervals are provided for 
#' all \eqn{x = 0, 1, 2 ..n} based on empirical Bayesian approach for Beta-Binomial model. 
#' Lower and Upper support values are needed to obtain the MLE of marginal likelihood 
#' for prior parameters.
#' @return A dataframe with 
#'  \item{x }{- Number of successes (positive samples)}
#'  \item{pomean }{ - Posterior mean}
#'  \item{LBAQ }{ - Lower limits of Quantile based intervals}
#'  \item{UBAQ }{ - Upper limits of Quantile based intervals}
#'  \item{LBAH }{ - Lower limits of HPD intervals}
#'  \item{UBAH }{ - Upper limits of HPD intervals}
#' @family Miscellaneous  functions for Bayesian method
#' @examples 
#' sL=runif(1,0,2)				#Lower and upper of Support for MLE optimization
#' sU=runif(1,sL,10)
#' n= 5; alp=0.05
#' empericalBA(n,alp,sL,sU) 
#' @references 
#' [1] 1998 Lehmann EL and Casella G
#' Theory of Point Estimation, 2nd ed Springer, New York
#' @export
##### DELTA_ALPHA, THETA,F-ERROR,POWER,FAILURE
empericalBA<-function(n,alp,sL,sU)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(sL)) stop("'Lower support' for MLE is missing")
  if (missing(sU)) stop("'Upper support' for MLE is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<0 ) stop("'n' has to be greater than or equal to 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(sL) != "integer") & (class(sL) != "numeric") || length(sL) >1 || sL<=0 ) stop("'sL' has to be a number greater than 0")
  if ((class(sU) != "integer") & (class(sU) != "numeric") || length(sU) >1 || sU<=sL ) stop("'sU' has to be a number greater than sL")
  
    x=0:n
    k=n+1
    ####INITIALIZATIONS
    av=0
    bv=0
    LEBAQ=0
    UEBAQ=0
    LEBAH=0
    UEBAH=0
    pomean=0
    #####Finding MLE of marginal likelihood for prior parameters in Binomial-Beta model
    for(i in 1:k)
    {
      likelhd = function(y) 
      {
        a<-y[1]
        b<-y[2]
        -(choose(n,x[i])/beta(a,b))*beta(x[i]+a,n-x[i]+b)
      }
      mle=optim(c(0.1,0.1), likelhd, NULL,method = "L-BFGS-B",lower = sL, upper = sU, hessian = TRUE)$par
      
      av[i]=mle[1]
      bv[i]=mle[2]
    }
   # library(TeachingDemos)				#To get HPDs
    
    #####Bayesian Posterior Estimates
    for(i in 1:k)
    {
      #Posterior mean
      pomean[i]=(x[i]+av[i])/(n+av[i]+bv[i])
      #Quantile Based Intervals
      LEBAQ[i]=qbeta(alp/2,x[i]+av[i],n-x[i]+bv[i])
      UEBAQ[i]=qbeta(1-(alp/2),x[i]+av[i],n-x[i]+bv[i])
      
      LEBAH[i]=TeachingDemos::hpd(qbeta,shape1=x[i]+av[i],shape2=n-x[i]+bv[i],conf=1-alp)[1]
      UEBAH[i]=TeachingDemos::hpd(qbeta,shape1=x[i]+av[i],shape2=n-x[i]+bv[i],conf=1-alp)[2]
    }
    return(round(data.frame(x,pomean,LEBAQ,UEBAQ,LEBAH,UEBAH),4))
  }
  