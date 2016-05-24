#'  Wald method of CI estimation
#' @param x - Number of successes
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @details  Wald-type interval that results from inverting
#' large-sample test and evaluates standard errors at maximum likelihood estimates
#' for the given \code{x} and \code{n}
#' @return A dataframe with
#'  \item{x}{  Number of successes (positive samples)}
#'  \item{LWDx }{   Wald Lower limit}
#'  \item{UWDx }{   Wald Upper Limit}
#'  \item{LABB }{   Wald Lower Abberation}
#'  \item{UABB }{   Wald Upper Abberation}
#'  \item{ZWI }{   Zero Width Interval}
#' @family Base methods of CI estimation given x & n
#' @seealso \code{\link{prop.test} and \link{binom.test}} for equivalent base Stats R functionality,
#'    \code{\link[binom]{binom.confint}}  provides similar functionality for 11 methods,
#'    \code{\link[PropCIs]{wald2ci}} which provides multiple functions for CI calculation ,
#'    \code{\link[BlakerCI]{binom.blaker.limits}} which calculates Blaker CI which is not covered here and
#'    \code{\link[prevalence]{propCI}} which provides similar functionality.
#' @examples
#' x=5; n=5; alp=0.05;
#' ciWDx(x,n,alp)
#' @references
#' [1] 1993 Vollset SE.
#' Confidence intervals for a binomial proportion.
#' Statistics in Medicine: 12; 809 - 824.
#'
#' [2] 1998 Agresti A and Coull BA.
#' Approximate is better than "Exact" for interval estimation of binomial proportions.
#' The American Statistician: 52; 119 - 126.
#'
#' [3] 1998 Newcombe RG.
#' Two-sided confidence intervals for the single proportion: Comparison of seven methods.
#' Statistics in Medicine: 17; 857 - 872.
#'
#' [4] 2001 Brown LD, Cai TT and DasGupta A.
#' Interval estimation for a binomial proportion.
#' Statistical Science: 16; 101 - 133.
#'
#' [5] 2002 Pan W.
#' Approximate confidence intervals for one proportion and difference of two proportions
#' Computational Statistics and Data Analysis 40, 128, 143-157.
#'
#' [6] 2008 Pires, A.M., Amado, C.
#' Interval Estimators for a Binomial Proportion: Comparison of Twenty Methods.
#' REVSTAT - Statistical Journal, 6, 165-197.
#'
#' [7] 2014 Martin Andres, A. and Alvarez Hernandez, M.
#' Two-tailed asymptotic inferences for a proportion.
#' Journal of Applied Statistics, 41, 7, 1516-1529
#' @export
#1.WALD
ciWDx<-function(x,n,alp)
{
  if (missing(x)) stop("'x' is missing")
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (((class(x) != "integer") & (class(x) != "numeric")) || (x<0) || x>n || length(x)>1) stop("'x' has to be a positive integer between 0 and n")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")


###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#WALD METHOD
pW=x/n
qW=1-(x/n)
seW=sqrt(pW*qW/n)
LWDx=pW-(cv*seW)
UWDx=pW+(cv*seW)
if(LWDx<0) LABB="YES" else LABB="NO"
if(LWDx<0) LWDx=0
if(UWDx>1) UABB="YES" else  UABB="NO"
if(UWDx>1) UWDx=1
if(UWDx-LWDx==0)ZWI="YES" else ZWI="NO"
return(data.frame(x,LWDx,UWDx,LABB,UABB,ZWI))
}
##############################################################################
#' Base Score method of CI estimation
#' @param x - Number of successes
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @details  A score test approach based on inverting the test
#' with standard error evaluated at the null hypothesis is due to
#' Wilson for the given \code{x} and \code{n}
#' @return A dataframe with
#'  \item{x}{  Number of successes (positive samples)}
#'  \item{LSCx }{   Score Lower limit}
#'  \item{USCx }{   Score Upper Limit}
#'  \item{LABB }{   Score Lower Abberation}
#'  \item{UABB }{   Score Upper Abberation}
#'  \item{ZWI }{   Zero Width Interval}
#' @family  Base methods of CI estimation given x & n
#' @seealso \code{\link{prop.test} and \link{binom.test}} for equivalent base Stats R functionality,
#'    \code{\link[binom]{binom.confint}}  provides similar functionality for 11 methods,
#'    \code{\link[PropCIs]{wald2ci}} which provides multiple functions for CI calculation ,
#'    \code{\link[BlakerCI]{binom.blaker.limits}} which calculates Blaker CI which is not covered here and
#'    \code{\link[prevalence]{propCI}} which provides similar functionality.
#' @examples
#' x=5; n=5; alp=0.05
#' ciSCx(x,n,alp)
#' @references
#' [1] 1993 Vollset SE.
#' Confidence intervals for a binomial proportion.
#' Statistics in Medicine: 12; 809 - 824.
#'
#' [2] 1998 Agresti A and Coull BA.
#' Approximate is better than "Exact" for interval estimation of binomial proportions.
#' The American Statistician: 52; 119 - 126.
#'
#' [3] 1998 Newcombe RG.
#' Two-sided confidence intervals for the single proportion: Comparison of seven methods.
#' Statistics in Medicine: 17; 857 - 872.
#'
#' [4] 2001 Brown LD, Cai TT and DasGupta A.
#' Interval estimation for a binomial proportion.
#' Statistical Science: 16; 101 - 133.
#'
#' [5] 2002 Pan W.
#' Approximate confidence intervals for one proportion and difference of two proportions
#' Computational Statistics and Data Analysis 40, 128, 143-157.
#'
#' [6] 2008 Pires, A.M., Amado, C.
#' Interval Estimators for a Binomial Proportion: Comparison of Twenty Methods.
#' REVSTAT - Statistical Journal, 6, 165-197.
#'
#' [7] 2014 Martin Andres, A. and Alvarez Hernandez, M.
#' Two-tailed asymptotic inferences for a proportion.
#' Journal of Applied Statistics, 41, 7, 1516-1529
#' @export
#2.SCORE
ciSCx<-function(x,n,alp)
{
  if (missing(x)) stop("'x' is missing")
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (((class(x) != "integer") & (class(x) != "numeric")) || (x<0) || x>n || length(x)>1) stop("'x' has to be a positive integer between 0 and n")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")

###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
cv1=(cv^2)/(2*n)
cv2=(cv/(2*n))^2

#SCORE (WILSON) METHOD
pS=x/n
qS=1-(x/n)
seS=sqrt((pS*qS/n)+cv2)
LSCx=(n/(n+(cv)^2))*((pS+cv1)-(cv*seS))
USCx=(n/(n+(cv)^2))*((pS+cv1)+(cv*seS))
if(LSCx<0) LABB="YES" else LABB="NO"
if(LSCx<0) LSCx=0

if(USCx>1) UABB="YES" else  UABB="NO"
if(USCx>1) USCx=1

if(USCx-LSCx==0)ZWI="YES" else ZWI="NO"

return(data.frame(x,LSCx,USCx,LABB,UABB,ZWI))
}
##############################################################################
#' Base ArcSine method of CI estimation
#' @param x - Number of successes
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @details  Wald-type interval for the given x and n using the arcsine
#' transformation of the parameter \code{p}; that is based on the normal approximation
#' for \eqn{sin^{-1}(p) }
#' @return A dataframe with
#'  \item{x}{  Number of successes (positive samples)}
#'  \item{LASx }{   ArcSine Lower limit}
#'  \item{UASx }{   ArcSine Upper Limit}
#'  \item{LABB }{   ArcSine Lower Abberation}
#'  \item{UABB }{   ArcSine Upper Abberation}
#'  \item{ZWI }{   Zero Width Interval}
#' @family Base methods of CI estimation given x & n
#' @seealso \code{\link{prop.test} and \link{binom.test}} for equivalent base Stats R functionality,
#'    \code{\link[binom]{binom.confint}}  provides similar functionality for 11 methods,
#'    \code{\link[PropCIs]{wald2ci}} which provides multiple functions for CI calculation ,
#'    \code{\link[BlakerCI]{binom.blaker.limits}} which calculates Blaker CI which is not covered here and
#'    \code{\link[prevalence]{propCI}} which provides similar functionality.
#' @examples
#' x=5; n=5; alp=0.05
#' ciASx(x,n,alp)
#' @references
#' [1] 1993 Vollset SE.
#' Confidence intervals for a binomial proportion.
#' Statistics in Medicine: 12; 809 - 824.
#'
#' [2] 1998 Agresti A and Coull BA.
#' Approximate is better than "Exact" for interval estimation of binomial proportions.
#' The American Statistician: 52; 119 - 126.
#'
#' [3] 1998 Newcombe RG.
#' Two-sided confidence intervals for the single proportion: Comparison of seven methods.
#' Statistics in Medicine: 17; 857 - 872.
#'
#' [4] 2001 Brown LD, Cai TT and DasGupta A.
#' Interval estimation for a binomial proportion.
#' Statistical Science: 16; 101 - 133.
#'
#' [5] 2002 Pan W.
#' Approximate confidence intervals for one proportion and difference of two proportions
#' Computational Statistics and Data Analysis 40, 128, 143-157.
#'
#' [6] 2008 Pires, A.M., Amado, C.
#' Interval Estimators for a Binomial Proportion: Comparison of Twenty Methods.
#' REVSTAT - Statistical Journal, 6, 165-197.
#'
#' [7] 2014 Martin Andres, A. and Alvarez Hernandez, M.
#' Two-tailed asymptotic inferences for a proportion.
#' Journal of Applied Statistics, 41, 7, 1516-1529
#' @export
#3.ARCSINE
ciASx<-function(x,n,alp)
{
  if (missing(x)) stop("'x' is missing")
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (((class(x) != "integer") & (class(x) != "numeric")) || (x<0) || x>n || length(x)>1) stop("'x' has to be a positive integer between 0 and n")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")

cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#ARC-SINE METHOD
pA=x/n
#qA=1-pA
seA=cv/sqrt(4*n)
LASx=(sin(asin(sqrt(pA))-seA))^2
UASx=(sin(asin(sqrt(pA))+seA))^2

if(LASx<0) LABB="YES" else LABB="NO"
if(LASx<0) LASx=0

if(UASx>1) UABB="YES" else  UABB="NO"
if(UASx>1) UASx=1

if(UASx-LASx==0)ZWI="YES" else ZWI="NO"
return(data.frame(x,LASx,UASx,LABB,UABB,ZWI))
}
##############################################################################
#' Likelihood Ratio method of CI estimation
#' @param x - Number of successes
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @details  Likelihood ratio limits for the given \code{x} and \code{n}
#' obtained as the solution to the equation in \code{p} formed as logarithm of
#' ratio between binomial likelihood at sample proportion and that of over
#' all possible parameters
#' @return A dataframe with
##'  \item{x }{- Number of successes (positive samples)}
##'  \item{LLRx }{ - Likelihood Ratio Lower limit}
##'  \item{ULRx }{ - Likelihood Ratio Upper Limit}
##'  \item{LABB }{ - Likelihood Ratio Lower Abberation}
##'  \item{UABB }{ - Likelihood Ratio Upper Abberation}
##'  \item{ZWI }{ - Zero Width Interval}
#' @family  Base methods of CI estimation given x & n
#' @seealso \code{\link{prop.test} and \link{binom.test}} for equivalent base Stats R functionality,
#'    \code{\link[binom]{binom.confint}}  provides similar functionality for 11 methods,
#'    \code{\link[PropCIs]{wald2ci}} which provides multiple functions for CI calculation ,
#'    \code{\link[BlakerCI]{binom.blaker.limits}} which calculates Blaker CI which is not covered here and
#'    \code{\link[prevalence]{propCI}} which provides similar functionality.
#' @examples
#' x=5; n=5; alp=0.05
#' ciLRx(x,n,alp)
#' @references
#' [1] 1993 Vollset SE.
#' Confidence intervals for a binomial proportion.
#' Statistics in Medicine: 12; 809 - 824.
#'
#' [2] 1998 Agresti A and Coull BA.
#' Approximate is better than "Exact" for interval estimation of binomial proportions.
#' The American Statistician: 52; 119 - 126.
#'
#' [3] 1998 Newcombe RG.
#' Two-sided confidence intervals for the single proportion: Comparison of seven methods.
#' Statistics in Medicine: 17; 857 - 872.
#'
#' [4] 2001 Brown LD, Cai TT and DasGupta A.
#' Interval estimation for a binomial proportion.
#' Statistical Science: 16; 101 - 133.
#'
#' [5] 2002 Pan W.
#' Approximate confidence intervals for one proportion and difference of two proportions
#' Computational Statistics and Data Analysis 40, 128, 143-157.
#'
#' [6] 2008 Pires, A.M., Amado, C.
#' Interval Estimators for a Binomial Proportion: Comparison of Twenty Methods.
#' REVSTAT - Statistical Journal, 6, 165-197.
#'
#' [7] 2014 Martin Andres, A. and Alvarez Hernandez, M.
#' Two-tailed asymptotic inferences for a proportion.
#' Journal of Applied Statistics, 41, 7, 1516-1529
#' @export
#4.LIKELIHOOD RATIO
ciLRx<-function(x,n,alp)
{
  if (missing(x)) stop("'x' is missing")
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (((class(x) != "integer") & (class(x) != "numeric")) || (x<0) || x>n || length(x)>1) stop("'x' has to be a positive integer between 0 and n")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")

####INPUT
y=x
###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#LIKELIHOOD-RATIO METHOD
likelhd = function(p) dbinom(y,n,p)
loglik = function(p) dbinom(y,n,p,log=TRUE)
mle=optimize(likelhd,c(0,1),maximum=TRUE)$maximum
cutoff=loglik(mle)-(cv^2/2)
loglik.optim=function(p){abs(cutoff-loglik(p))}
LLRx=optimize(loglik.optim, c(0,mle))$minimum
ULRx=optimize(loglik.optim, c(mle,1))$minimum

if(LLRx<0) LABB="YES" else LABB="NO"
if(LLRx<0) LLRx=0

if(ULRx>1) UABB="YES" else  UABB="NO"
if(ULRx>1) ULRx=1

if(ULRx-LLRx==0)ZWI="YES" else ZWI="NO"
return(data.frame(x,LLRx,ULRx,LABB,UABB,ZWI))
}
##############################################################################
#' Exact method of CI estimation
#' @param n - Number of trials
#' @param x - Number of sucess
#' @param alp - Alpha value (significance level required)
#' @param e - Exact method indicator  in [0, 1] {1: Clopper Pearson, 0.5: Mid P}
#' @details  Confidence interval for \code{p} (for the given \code{x} and \code{n}),
#' based on inverting equal-tailed binomial tests with null hypothesis \eqn{H0: p = p0 }
#' and calculated from the cumulative binomial distribution.
#' Exact two sided P-value is usually calculated as \deqn{P= 2[e*Pr(X = x) + min({Pr(X < x), Pr(X > x)})]}
#'  where probabilities are found at null value of \code{p} and \eqn{0 <= e <= 1.}
#' @return A dataframe with
##'  \item{x }{- Number of successes (positive samples)}
##'  \item{LEXx }{ - Exact Lower limit}
##'  \item{UEXx }{ - Exact Upper Limit}
#' @family  Base methods of CI estimation given x & n
#' @seealso \code{\link{prop.test} and \link{binom.test}} for equivalent base Stats R functionality,
#'    \code{\link[binom]{binom.confint}}  provides similar functionality for 11 methods,
#'    \code{\link[PropCIs]{wald2ci}} which provides multiple functions for CI calculation ,
#'    \code{\link[BlakerCI]{binom.blaker.limits}} which calculates Blaker CI which is not covered here and
#'    \code{\link[prevalence]{propCI}} which provides similar functionality.
#' @examples
#' x=5; n=5; alp=0.05;e=0.5
#' ciEXx(x,n,alp,e) #Mid-p
#' x=5; n=5; alp=0.05;e=1 #Clopper Pearson
#' ciEXx(x,n,alp,e)
#' x=5; n=5; alp=0.05;e=c(0.1,0.5,0.95,1) #Range including Mid-p and Clopper-Pearson
#' ciEXx(x,n,alp,e)
#' @references
#' [1] 1993 Vollset SE.
#' Confidence intervals for a binomial proportion.
#' Statistics in Medicine: 12; 809 - 824.
#'
#' [2] 1998 Agresti A and Coull BA.
#' Approximate is better than "Exact" for interval estimation of binomial proportions.
#' The American Statistician: 52; 119 - 126.
#'
#' [3] 1998 Newcombe RG.
#' Two-sided confidence intervals for the single proportion: Comparison of seven methods.
#' Statistics in Medicine: 17; 857 - 872.
#'
#' [4] 2001 Brown LD, Cai TT and DasGupta A.
#' Interval estimation for a binomial proportion.
#' Statistical Science: 16; 101 - 133.
#'
#' [5] 2002 Pan W.
#' Approximate confidence intervals for one proportion and difference of two proportions
#' Computational Statistics and Data Analysis 40, 128, 143-157.
#'
#' [6] 2008 Pires, A.M., Amado, C.
#' Interval Estimators for a Binomial Proportion: Comparison of Twenty Methods.
#' REVSTAT - Statistical Journal, 6, 165-197.
#'
#' [7] 2014 Martin Andres, A. and Alvarez Hernandez, M.
#' Two-tailed asymptotic inferences for a proportion.
#' Journal of Applied Statistics, 41, 7, 1516-1529
#' @export
#5.EXACT METHOD
ciEXx=function(x,n,alp,e)
{
  if (missing(x)) stop("'x' is missing")
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(e)) stop("'e' is missing")
  if (((class(x) != "integer") & (class(x) != "numeric")) || (x<0) || x>n || length(x)>1) stop("'x' has to be a positive integer between 0 and n")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(e) != "integer") & (class(e) != "numeric") || any(e>1) || any(e<0)) stop("'e' has to be between 0 and 1")
  if (length(e)>10) stop("'e' can have only 10 intervals")

  nvar=length(e)

  res <- data.frame()

  for(i in 1:nvar)
  {
    lu=lufn103(x,n,alp,e[i])
    res <- rbind(res,lu)
  }
  return(res)
}
lufn103<-function(x,n,alp,e)
{
  LEXx=0
  UEXx=0
  LABB=0
  UABB=0
  ZWI=0

  LEXx=exlim103l(x,n,alp,e)
  UEXx=exlim103u(x,n,alp,e)
  if(LEXx<0) LABB="YES" else LABB="NO"
  if(LEXx<0) LEXx=0

  if(UEXx>1) UABB="YES" else UABB="NO"
  if(UEXx>1) UEXx=1

  if(UEXx-LEXx==0) ZWI="YES" else ZWI="NO"

    resx=data.frame(x,LEXx,UEXx,LABB, UABB, ZWI, e)
    return(resx)
}
exlim103l=function(x,n,alp,e)
{
  if(x==0)
  {
    LEX = 0
  } else if(x==n){
    LEX= (alp/(2*e))^(1/n)
  }else
  {
    z=x-1
    y=0:z
    f1=function(p) (1-e)*dbinom(x,n,p)+sum(dbinom(y,n,p))-(1-(alp/2))
    LEX= uniroot(f1,c(0,1))$root
  }
  return(LEX)
}
exlim103u=function(x,n,alp,e)
{
  if(x==0)
  {
    UEX = 1-((alp/(2*e))^(1/n))
  } else if(x==n){
    UEX = 1
  }else
  {
    z=x-1
    y=0:z
    f2= function(p) e*dbinom(x,n,p)+sum(dbinom(y,n,p))-(alp/2)
    UEX =uniroot(f2,c(0,1))$root
  }
  return(UEX)
}
######################################################################
#' Wald-T method of CI estimation
#' @param x - Number of successes
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @details  An approximate method based on a t_approximation of the
#' standardized point estimator for the given \code{x} and \code{n}; that is the point estimator
#' divided by its estimated standard error. Essential boundary modification is when \eqn{x = 0 or n},
#' \deqn{\hat{p}=\frac{(x+2)}{(n+4)}{phat=(x+2)/(n+4)}}
#' @return A dataframe with
#'  \item{x}{  Number of successes (positive samples)}
#'  \item{LTWx }{   Wald-T Lower limit}
#'  \item{UTWx }{   Wald-T Upper Limit}
#'  \item{LABB }{   Wald-T Lower Abberation}
#'  \item{UABB }{   Wald-T Upper Abberation}
#'  \item{ZWI }{   Zero Width Interval}
#' @family Base methods of CI estimation given x & n
#' @seealso \code{\link{prop.test} and \link{binom.test}} for equivalent base Stats R functionality,
#'    \code{\link[binom]{binom.confint}}  provides similar functionality for 11 methods,
#'    \code{\link[PropCIs]{wald2ci}} which provides multiple functions for CI calculation ,
#'    \code{\link[BlakerCI]{binom.blaker.limits}} which calculates Blaker CI which is not covered here and
#'    \code{\link[prevalence]{propCI}} which provides similar functionality.
#' @examples
#' x=5; n=5; alp=0.05
#' ciTWx(x,n,alp)
#' @references
#' [1] 1993 Vollset SE.
#' Confidence intervals for a binomial proportion.
#' Statistics in Medicine: 12; 809 - 824.
#'
#' [2] 1998 Agresti A and Coull BA.
#' Approximate is better than "Exact" for interval estimation of binomial proportions.
#' The American Statistician: 52; 119 - 126.
#'
#' [3] 1998 Newcombe RG.
#' Two-sided confidence intervals for the single proportion: Comparison of seven methods.
#' Statistics in Medicine: 17; 857 - 872.
#'
#' [4] 2001 Brown LD, Cai TT and DasGupta A.
#' Interval estimation for a binomial proportion.
#' Statistical Science: 16; 101 - 133.
#'
#' [5] 2002 Pan W.
#' Approximate confidence intervals for one proportion and difference of two proportions
#' Computational Statistics and Data Analysis 40, 128, 143-157.
#'
#' [6] 2008 Pires, A.M., Amado, C.
#' Interval Estimators for a Binomial Proportion: Comparison of Twenty Methods.
#' REVSTAT - Statistical Journal, 6, 165-197.
#'
#' [7] 2014 Martin Andres, A. and Alvarez Hernandez, M.
#' Two-tailed asymptotic inferences for a proportion.
#' Journal of Applied Statistics, 41, 7, 1516-1529
#' @export
#6.WALD-T
ciTWx<-function(x,n,alp)
{
  if (missing(x)) stop("'x' is missing")
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (((class(x) != "integer") & (class(x) != "numeric")) || (x<0) || x>n || length(x)>1) stop("'x' has to be a positive integer between 0 and n")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")

#MODIFIED_t-WALD METHOD
if(x==0||x==n)
{
pTWx=(x+2)/(n+4)
#qTWx=1-pTWx
}else
{
pTWx=x/n
#qTWx=1-pTWx
}
f1=function(p,n) p*(1-p)/n
f2=function(p,n) (p*(1-p)/(n^3))+(p+((6*n)-7)*(p^2)+(4*(n-1)*(n-3)*(p^3))-(2*(n-1)*((2*n)-3)*(p^4)))/(n^5)-(2*(p+((2*n)-3)*(p^2)-2*(n-1)*(p^3)))/(n^4)
DOFx=2*((f1(pTWx,n))^2)/f2(pTWx,n)
cvx=qt(1-(alp/2), df=DOFx)
seTWx=cvx*sqrt(f1(pTWx,n))
LTWx=pTWx-(seTWx)
UTWx=pTWx+(seTWx)

if(LTWx<0) LABB="YES" else LABB="NO"
if(LTWx<0) LTWx=0

if(UTWx>1) UABB="YES" else  UABB="NO"
if(UTWx>1) UTWx=1

if(UTWx-LTWx==0)ZWI="YES" else ZWI="NO"
return(data.frame(x,LTWx,UTWx,LABB,UABB,ZWI))
}
#####################################################################
#' Logit Wald method of CI estimation
#' @param x - Number of successes
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @details  Wald-type interval for all the given \code{x} and \code{n}
#' based on the logit transformation of \code{p};
#' that is normal approximation for \eqn{log(p/1-p)}
#' @return A dataframe with
#'  \item{x }{- Number of successes (positive samples)}
#'  \item{LLTx }{ - Logit Wald Lower limit}
#'  \item{ULTx }{ - Logit Wald Upper Limit}
#'  \item{LABB }{ - Logit Wald Lower Abberation}
#'  \item{UABB }{ - Logit Wald Upper Abberation}
#'  \item{ZWI }{ - Zero Width Interval}
#' @family Base methods of CI estimation given x & n
#' @seealso \code{\link{prop.test} and \link{binom.test}} for equivalent base Stats R functionality,
#'    \code{\link[binom]{binom.confint}}  provides similar functionality for 11 methods,
#'    \code{\link[PropCIs]{wald2ci}} which provides multiple functions for CI calculation ,
#'    \code{\link[BlakerCI]{binom.blaker.limits}} which calculates Blaker CI which is not covered here and
#'    \code{\link[prevalence]{propCI}} which provides similar functionality.
#' @examples
#' x=5; n=5; alp=0.05
#' ciLTx(x,n,alp)
#' @references
#' [1] 1993 Vollset SE.
#' Confidence intervals for a binomial proportion.
#' Statistics in Medicine: 12; 809 - 824.
#'
#' [2] 1998 Agresti A and Coull BA.
#' Approximate is better than "Exact" for interval estimation of binomial proportions.
#' The American Statistician: 52; 119 - 126.
#'
#' [3] 1998 Newcombe RG.
#' Two-sided confidence intervals for the single proportion: Comparison of seven methods.
#' Statistics in Medicine: 17; 857 - 872.
#'
#' [4] 2001 Brown LD, Cai TT and DasGupta A.
#' Interval estimation for a binomial proportion.
#' Statistical Science: 16; 101 - 133.
#'
#' [5] 2002 Pan W.
#' Approximate confidence intervals for one proportion and difference of two proportions
#' Computational Statistics and Data Analysis 40, 128, 143-157.
#'
#' [6] 2008 Pires, A.M., Amado, C.
#' Interval Estimators for a Binomial Proportion: Comparison of Twenty Methods.
#' REVSTAT - Statistical Journal, 6, 165-197.
#'
#' [7] 2014 Martin Andres, A. and Alvarez Hernandez, M.
#' Two-tailed asymptotic inferences for a proportion.
#' Journal of Applied Statistics, 41, 7, 1516-1529
#' @export
#7.LOGIT-WALD
ciLTx<-function(x,n,alp)
{
  if (missing(x)) stop("'x' is missing")
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (((class(x) != "integer") & (class(x) != "numeric")) || (x<0) || x>n || length(x)>1) stop("'x' has to be a positive integer between 0 and n")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")

####INITIALIZATIONS
pLTx=0
qLTx=0
seLTx=0
lgitx=0
LLTx=0
ULTx=0
LABB=0
UABB=0
ZWI=0
###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#LOGIT-WALD METHOD
if(x==0)
{
pLTx=0
qLTx=1
LLTx = 0
ULTx = 1-((alp/2)^(1/n))
}
if(x==n)
{
pLTx=1
qLTx=0
LLTx= (alp/2)^(1/n)
ULTx=1
}
if(0<x&&x<n)
{
pLTx=x/n
qLTx=1-pLTx
lgitx=log(pLTx/qLTx)
seLTx=sqrt(pLTx*qLTx*n)
LLTx=1/(1+exp(-lgitx+(cv/seLTx)))
ULTx=1/(1+exp(-lgitx-(cv/seLTx)))
}

if(LLTx<0) LABB="YES" else LABB="NO"
if(LLTx<0) LLTx=0

if(ULTx>1) UABB="YES" else UABB="NO"
if(ULTx>1) ULTx=1

if(ULTx-LLTx==0)ZWI="YES" else ZWI="NO"

return(data.frame(x,LLTx,ULTx,LABB,UABB,ZWI))
}
#####################################################################################
#' Bayesian method of CI estimation  with Beta prior distribution
#' @param x - Number of sucess
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Shape parameter 1 for prior Beta distribution in Bayesian model. Can also be a vector of length n+1 priors.
#' @param b - Shape parameter 2 for prior Beta distribution in Bayesian model. Can also be a vector of length n+1 priors.
#' @details  Highest Probability Density (HPD) and two tailed intervals are
#' provided for the given \code{x} and \code{n}. based on the conjugate prior \eqn{\beta(a, b)}
#' for the probability of success \code{p} of the binomial distribution so that the
#' posterior is \eqn{\beta(x + a, n - x + b)}.
#' @return A dataframe with
##'  \item{x }{- Number of successes (positive samples)}
##'  \item{LBAQx }{ - Lower limits of Quantile based intervals}
##'  \item{UBAQx }{ - Upper limits of Quantile based intervals}
##'  \item{LBAHx }{ - Lower limits of HPD intervals}
##'  \item{UBAHx }{ - Upper limits of HPD intervals}
#' @family Base methods of CI estimation given x & n
#' @seealso \code{\link{prop.test} and \link{binom.test}} for equivalent base Stats R functionality,
#'    \code{\link[binom]{binom.confint}}  provides similar functionality for 11 methods,
#'    \code{\link[PropCIs]{wald2ci}} which provides multiple functions for CI calculation ,
#'    \code{\link[BlakerCI]{binom.blaker.limits}} which calculates Blaker CI which is not covered here and
#'    \code{\link[prevalence]{propCI}} which provides similar functionality.
#' @examples
#' x=5; n=5; alp=0.05; a=0.5;b=0.5;
#' ciBAx(x,n,alp,a,b)
#' x= 5; n=5; alp=0.05; a=c(0.5,2,1,1,2,0.5);b=c(0.5,2,1,1,2,0.5)
#' ciBAx(x,n,alp,a,b)
#' @references
#' [1] 1993 Vollset SE.
#' Confidence intervals for a binomial proportion.
#' Statistics in Medicine: 12; 809 - 824.
#'
#' [2] 1998 Agresti A and Coull BA.
#' Approximate is better than "Exact" for interval estimation of binomial proportions.
#' The American Statistician: 52; 119 - 126.
#'
#' [3] 1998 Newcombe RG.
#' Two-sided confidence intervals for the single proportion: Comparison of seven methods.
#' Statistics in Medicine: 17; 857 - 872.
#'
#' [4] 2001 Brown LD, Cai TT and DasGupta A.
#' Interval estimation for a binomial proportion.
#' Statistical Science: 16; 101 - 133.
#'
#' [5] 2002 Pan W.
#' Approximate confidence intervals for one proportion and difference of two proportions
#' Computational Statistics and Data Analysis 40, 128, 143-157.
#'
#' [6] 2008 Pires, A.M., Amado, C.
#' Interval Estimators for a Binomial Proportion: Comparison of Twenty Methods.
#' REVSTAT - Statistical Journal, 6, 165-197.
#'
#' [7] 2014 Martin Andres, A. and Alvarez Hernandez, M.
#' Two-tailed asymptotic inferences for a proportion.
#' Journal of Applied Statistics, 41, 7, 1516-1529
#' @export
#8. BAYESIAN
ciBAx<-function(x,n,alp,a,b)
{
  if (missing(x)) stop("'x' is missing")
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (((class(x) != "integer") & (class(x) != "numeric")) || (x<0) || x>n || length(x)>1) stop("'x' has to be a positive integer between 0 and n")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")

  if(length(a)>1 || length(b)>1) {
    if ((class(a) != "integer") & (class(a) != "numeric") || any(a < 0)) stop("'a' has to be a set of positive numeric vectors")
    if ((class(b) != "integer") & (class(b) != "numeric") || any(b < 0)) stop("'b' has to be a set of positive numeric vectors")
    if (length(a) <=  n || length(b) <=  n ) stop("'a' and 'b' vectors have to be equal to or greater than length n+1")
    BAdfx=ciBADx(x,n,alp,a,b)
    }
  else
    {
  if (a<0 ) stop("'a' has to be greater than or equal to 0")
  if (b<0 ) stop("'b' has to be greater than or equal to 0")

##############
#library(TeachingDemos)				#To get HPDs
#Quantile Based Intervals
LBAQx=qbeta(alp/2,x+a,n-x+b)
UBAQx=qbeta(1-(alp/2),x+a,n-x+b)

LBAHx=TeachingDemos::hpd(qbeta,shape1=x+a,shape2=n-x+b,conf=1-alp)[1]
UBAHx=TeachingDemos::hpd(qbeta,shape1=x+a,shape2=n-x+b,conf=1-alp)[2]

BAdfx=data.frame(x,LBAQx,UBAQx,LBAHx,UBAHx)
}

return(BAdfx)
}
#######################################################################
ciBADx<-function(x,n,alp,a,b)
{

  k=n+1
  ####INITIALIZATIONS
  LBAQx=0
  UBAQx=0
  LBAHx=0
  UBAHx=0
  ##############

  for(i in 1:k)
  {
    #Quantile Based Intervals
    LBAQx[i]=qbeta(alp/2,x+a[i],n-x+b[i])
    UBAQx[i]=qbeta(1-(alp/2),x+a[i],n-x+b[i])

    LBAHx[i]=TeachingDemos::hpd(qbeta,shape1=x+a[i],shape2=n-x+b[i],conf=1-alp)[1]
    UBAHx[i]=TeachingDemos::hpd(qbeta,shape1=x+a[i],shape2=n-x+b[i],conf=1-alp)[2]

  }

  return(data.frame(x,LBAQx,UBAQx,LBAHx,UBAHx))
}
#####################################################################################
#' Specific CI estimation of 6 base methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)
#' @param x - Number of sucess
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @details  The Confidence Interval of using 6 base methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)  for \code{n} given \code{alp} and \code{x}
#' @return A dataframe with
##'  \item{name }{- Name of the method}
##'  \item{x }{- Number of successes (positive samples)}
##'  \item{LLT }{ - Lower limit}
##'  \item{ULT }{ - Upper Limit}
##'  \item{LABB }{ - Lower Abberation}
##'  \item{UABB }{ - Upper Abberation}
##'  \item{ZWI }{ - Zero Width Interval}
#' @family Base methods of CI estimation given x & n
#' @seealso \code{\link{prop.test} and \link{binom.test}} for equivalent base Stats R functionality,
#'    \code{\link[binom]{binom.confint}}  provides similar functionality for 11 methods,
#'    \code{\link[PropCIs]{wald2ci}} which provides multiple functions for CI calculation ,
#'    \code{\link[BlakerCI]{binom.blaker.limits}} which calculates Blaker CI which is not covered here and
#'    \code{\link[prevalence]{propCI}} which provides similar functionality.
#' @examples
#' x= 5; n=5; alp=0.05;
#' ciAllx(x,n,alp)
#' @references
#' [1] 1993 Vollset SE.
#' Confidence intervals for a binomial proportion.
#' Statistics in Medicine: 12; 809 - 824.
#'
#' [2] 1998 Agresti A and Coull BA.
#' Approximate is better than "Exact" for interval estimation of binomial proportions.
#' The American Statistician: 52; 119 - 126.
#'
#' [3] 1998 Newcombe RG.
#' Two-sided confidence intervals for the single proportion: Comparison of seven methods.
#' Statistics in Medicine: 17; 857 - 872.
#'
#' [4] 2001 Brown LD, Cai TT and DasGupta A.
#' Interval estimation for a binomial proportion.
#' Statistical Science: 16; 101 - 133.
#'
#' [5] 2002 Pan W.
#' Approximate confidence intervals for one proportion and difference of two proportions
#' Computational Statistics and Data Analysis 40, 128, 143-157.
#'
#' [6] 2008 Pires, A.M., Amado, C.
#' Interval Estimators for a Binomial Proportion: Comparison of Twenty Methods.
#' REVSTAT - Statistical Journal, 6, 165-197.
#'
#' [7] 2014 Martin Andres, A. and Alvarez Hernandez, M.
#' Two-tailed asymptotic inferences for a proportion.
#' Journal of Applied Statistics, 41, 7, 1516-1529
#' @export
#9.All methods
ciAllx<-function(x,n,alp)
{
  if (missing(x)) stop("'x' is missing")
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (((class(x) != "integer") & (class(x) != "numeric")) || (x<0) || x>n || length(x)>1) stop("'x' has to be a positive integer between 0 and n")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")

  #### Calling functions and creating df
  WaldCI.df    = ciWDx(x,n,alp)
  ArcSineCI.df = ciASx(x,n,alp)
  LRCI.df      = ciLRx(x,n,alp)
  ScoreCI.df   = ciSCx(x,n,alp)
  WaldLCI.df   = ciLTx(x,n,alp)
  AdWaldCI.df  = ciTWx(x,n,alp)

  WaldCI.df$method    = as.factor("Wald")
  ArcSineCI.df$method = as.factor("ArcSine")
  LRCI.df$method      = as.factor("Likelihood")
  WaldLCI.df$method    = as.factor("Logit-Wald")
  ScoreCI.df$method   = as.factor("Score")
  AdWaldCI.df$method  = as.factor("Wald-T")

  Generic.1 = data.frame(method = WaldCI.df$method, x=WaldCI.df$x, LowerLimit = WaldCI.df$LWDx, UpperLimit = WaldCI.df$UWDx, LowerAbb = WaldCI.df$LABB, UpperAbb = WaldCI.df$UABB, ZWI = WaldCI.df$ZWI)
  Generic.2 = data.frame(method = ArcSineCI.df$method, x=ArcSineCI.df$x, LowerLimit = ArcSineCI.df$LASx, UpperLimit = ArcSineCI.df$UASx, LowerAbb = ArcSineCI.df$LABB, UpperAbb = ArcSineCI.df$UABB, ZWI = ArcSineCI.df$ZWI)
  Generic.3 = data.frame(method = LRCI.df$method, x=LRCI.df$x, LowerLimit = LRCI.df$LLRx, UpperLimit = LRCI.df$ULRx, LowerAbb = LRCI.df$LABB, UpperAbb = LRCI.df$UABB, ZWI = LRCI.df$ZWI)
  Generic.4 = data.frame(method = ScoreCI.df$method, x=ScoreCI.df$x, LowerLimit = ScoreCI.df$LSCx, UpperLimit = ScoreCI.df$USCx, LowerAbb = ScoreCI.df$LABB, UpperAbb = ScoreCI.df$UABB, ZWI = ScoreCI.df$ZWI)
  Generic.5 = data.frame(method = WaldLCI.df$method, x=WaldLCI.df$x, LowerLimit = WaldLCI.df$LLTx, UpperLimit = WaldLCI.df$ULTx, LowerAbb = WaldLCI.df$LABB, UpperAbb = WaldLCI.df$UABB, ZWI = WaldLCI.df$ZWI)
  Generic.6 = data.frame(method = AdWaldCI.df$method, x=AdWaldCI.df$x, LowerLimit = AdWaldCI.df$LTWx, UpperLimit = AdWaldCI.df$UTWx, LowerAbb = AdWaldCI.df$LABB, UpperAbb = AdWaldCI.df$UABB, ZWI = AdWaldCI.df$ZWI)

  Final.df= rbind(Generic.1,Generic.2,Generic.3,Generic.4,Generic.5, Generic.6)

  return(Final.df)
}
