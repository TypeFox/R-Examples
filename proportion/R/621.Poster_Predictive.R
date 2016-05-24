#' The Predicted probability - Bayesian approach
#' @param n - Number of trials from data
#' @param m - Future :Number of trials
#' @param a1 - Beta Prior Parameters for Bayesian estimation
#' @param a2 - Beta Prior Parameters for Bayesian estimation
#' @details  Computes posterior predictive probabilities for the required size of number of
#' trials \code{m} from the given number of trials \code{n} for the given parameters for Beta prior
#' distribution
#' @return A matrix of probability values between [0,1]
#'  \item{predicted_probability }{- The predicted probability}
#'  \item{0:n}{The number of columns based on the value of n}
#' @family Miscellaneous  functions for Bayesian method
#' @examples
#' n=10; m=5; a1=0.5; a2=0.5
#' probPRE(n,m,a1,a2)
#' @references
#' [1] 2002 Gelman A, Carlin  JB, Stern HS and Dunson DB
#' Bayesian Data Analysis, Chapman & Hall/CRC
#' @export
probPRE<-function(n,m,a1,a2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(m)) stop("'m' is missing")
  if (missing(a1)) stop("'a1' is missing")
  if (missing(a2)) stop("'a2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<0 ) stop("'n' has to be greater or equal to 0")
  if ((class(m) != "integer") & (class(m) != "numeric") || length(m) >1|| m<0 ) stop("'m' has to be greater or equal to 0")
  if ((class(a1) != "integer") & (class(a2) != "numeric") || length(a1) >1|| a1<=0 ) stop("'a1' has to be greater or equal to 0")
  if ((class(a2) != "integer") & (class(a2) != "numeric") || length(a2) >1|| a2<=0 ) stop("'a2' has to be greater or equal to 0")

  x=0:n
  xnew=0:m
  k1=n+1
  k2=m+1
  prepro=matrix(0,k2,k1)		#Predictive Probabilities

  for(j in 1:k2)
  {
    for(i in 1:k1)
    {
      prepro[j,i]=(choose(m,xnew[j]))/(beta(x[i]+a1,n-x[i]+a2))*beta(xnew[j]+x[i]+a1,m+n-xnew[j]-x[i]+a2)
    }
  }
  qq=matrix(c(xnew,prepro),m+1,n+2)
  colnames(qq)=c("xnew",0:n)

  return(predicted_probability=qq)
}
