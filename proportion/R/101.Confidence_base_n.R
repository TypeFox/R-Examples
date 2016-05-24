#'  Wald method of CI estimation
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @details  Wald-type interval that results from inverting large-sample test and evaluates standard errors at maximum likelihood estimates for all x = 0, 1, 2 ..n.
#' Calculates the Confidence Interval of \code{n} given \code{alp} along with lower and upper abberation.
#' @return A dataframe with
#'  \item{x}{  Number of successes (positive samples)}
#'  \item{LWD }{   Wald Lower limit}
#'  \item{UWD }{   Wald Upper Limit}
#'  \item{LABB }{   Wald Lower Abberation}
#'  \item{UABB }{   Wald Upper Abberation}
#'  \item{ZWI }{   Zero Width Interval}
#' @family Basic methods of CI estimation
#' @seealso \code{\link{prop.test} and \link{binom.test}} for equivalent base Stats R functionality,
#'    \code{\link[binom]{binom.confint}}  provides similar functionality for 11 methods,
#'    \code{\link[PropCIs]{wald2ci}} which provides multiple functions for CI calculation ,
#'    \code{\link[BlakerCI]{binom.blaker.limits}} which calculates Blaker CI which is not covered here and
#'    \code{\link[prevalence]{propCI}} which provides similar functionality.
#' @examples
#' n=5; alp=0.05
#' ciWD(n,alp)
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
ciWD<-function(n,alp)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")

  ####INPUT n
  x=0:n
  k=n+1
  ####INITIALIZATIONS
  pW=0
  qW=0
  seW=0
  LWD=0
  UWD=0
  LABB=0
  UABB=0
  ZWI=0
  ###CRITICAL VALUES
  cv=qnorm(1-(alp/2), mean = 0, sd = 1)
  #WALD METHOD
  for(i in 1:k)
  {
    pW[i]=x[i]/n
    qW[i]=1-(x[i]/n)
    seW[i]=sqrt(pW[i]*qW[i]/n)
    LWD[i]=pW[i]-(cv*seW[i])
    UWD[i]=pW[i]+(cv*seW[i])

    if(LWD[i]<0) LABB[i]="YES" else LABB[i]="NO"
    if(LWD[i]<0) LWD[i]=0

    if(UWD[i]>1) UABB[i]="YES" else UABB[i]="NO"
    if(UWD[i]>1) UWD[i]=1

    if(UWD[i]-LWD[i]==0)ZWI[i]="YES" else ZWI[i]="NO"
  }
  return(data.frame(x,LWD,UWD,LABB,UABB,ZWI))
}
##############################################################################
#'  Score method of CI estimation
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @details  A score test approach based on inverting the test with standard error evaluated at the null hypothesis is due to Wilson for all x = 0, 1, 2 ..n.
#  Calculates the Confidence Interval of \code{n} given \code{alp} along with lower and upper abberation.
#' @return A dataframe with
#'  \item{x }{- Number of successes (positive samples)}
#'  \item{LSC }{ - Score Lower limit}
#'  \item{USC }{ - Score Upper Limit}
#'  \item{LABB }{ - Score Lower Abberation}
#'  \item{UABB }{ - Score Upper Abberation}
#'  \item{ZWI }{ - Zero Width Interval}
#' @family Basic methods of CI estimation
#' @seealso \code{\link{prop.test} and \link{binom.test}} for equivalent base Stats R functionality,
#'    \code{\link[binom]{binom.confint}}  provides similar functionality for 11 methods,
#'    \code{\link[PropCIs]{wald2ci}} which provides multiple functions for CI calculation ,
#'    \code{\link[BlakerCI]{binom.blaker.limits}} which calculates Blaker CI which is not covered here and
#'    \code{\link[prevalence]{propCI}} which provides similar functionality.
#' @examples
#' n=5; alp=0.05
#' ciSC(n,alp)
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
ciSC<-function(n,alp)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")

  ####INPUT n
  x=0:n
  k=n+1
  ####INITIALIZATIONS
  pS=0
  qS=0
  seS=0
  LSC=0
  USC=0
  LABB=0
  UABB=0
  ZWI=0
  ###CRITICAL VALUES
  cv=qnorm(1-(alp/2), mean = 0, sd = 1)
  cv1=(cv^2)/(2*n)
  cv2=(cv/(2*n))^2

  #SCORE (WILSON) METHOD
  for(i in 1:k)
  {
    pS[i]=x[i]/n
    qS[i]=1-(x[i]/n)
    seS[i]=sqrt((pS[i]*qS[i]/n)+cv2)
    LSC[i]=(n/(n+(cv)^2))*((pS[i]+cv1)-(cv*seS[i]))
    USC[i]=(n/(n+(cv)^2))*((pS[i]+cv1)+(cv*seS[i]))

    if(LSC[i]<0) LABB[i]="YES" else LABB[i]="NO"
    if(LSC[i]<0) LSC[i]=0

    if(USC[i]>1) UABB[i]="YES" else UABB[i]="NO"
    if(USC[i]>1) USC[i]=1

    if(USC[i]-LSC[i]==0)ZWI[i]="YES" else ZWI[i]="NO"
  }
  return(data.frame(x,LSC,USC,LABB,UABB,ZWI))
}
##############################################################################
#'  ArcSine method of CI estimation
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @details  Wald-type interval for all x = 0, 1, 2 ..n using the arcsine transformation of the parameter p; that is based on the normal approximation for sin-1(p).
#' Calculates the Confidence Interval of \code{n} given \code{alp} along with lower and upper abberation.
#' @return A dataframe with
##'  \item{x }{- Number of successes (positive samples)}
##'  \item{LAS }{ - ArcSine Lower limit}
##'  \item{UAS }{ - ArcSine Upper Limit}
##'  \item{LABB }{ - ArcSine Lower Abberation}
##'  \item{UABB }{ - ArcSine Upper Abberation}
##'  \item{ZWI }{ - Zero Width Interval}
#' @family Basic methods of CI estimation
#' @seealso \code{\link{prop.test} and \link{binom.test}} for equivalent base Stats R functionality,
#'    \code{\link[binom]{binom.confint}}  provides similar functionality for 11 methods,
#'    \code{\link[PropCIs]{wald2ci}} which provides multiple functions for CI calculation ,
#'    \code{\link[BlakerCI]{binom.blaker.limits}} which calculates Blaker CI which is not covered here and
#'    \code{\link[prevalence]{propCI}} which provides similar functionality.
#' @examples
#' n=5; alp=0.05
#' ciAS(n,alp)
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
ciAS<-function(n,alp)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")

  ####INPUT
  x=0:n
  k=n+1
  ####INITIALIZATIONS
  pA=0
  qA=0
  seA=0
  LAS=0
  UAS=0
  LABB=0
  UABB=0
  ZWI=0

  cv=qnorm(1-(alp/2), mean = 0, sd = 1)
  #ARC-SINE METHOD
  for(i in 1:k)
  {
    pA[i]=x[i]/n
    qA[i]=1-pA[i]
    seA[i]=cv/sqrt(4*n)
    LAS[i]=(sin(asin(sqrt(pA[i]))-seA[i]))^2
    UAS[i]=(sin(asin(sqrt(pA[i]))+seA[i]))^2

    if(LAS[i]<0) LABB[i]="YES" else LABB[i]="NO"
    if(LAS[i]<0) LAS[i]=0

    if(UAS[i]>1) UABB[i]="YES" else UABB[i]="NO"
    if(UAS[i]>1) UAS[i]=1

    if(UAS[i]-LAS[i]==0)ZWI[i]="YES" else ZWI[i]="NO"
  }
  return(data.frame(x,LAS,UAS,LABB,UABB,ZWI))
}
##############################################################################

#'  Likelihood Ratio method of CI estimation
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @details  Likelihood ratio limits for all x = 0, 1, 2 ..n obtained as the solution to the
#' equation in p formed as logarithm of ratio between binomial likelihood at sample proportion
#' and that of over all possible parameters.
#' Calculates the Confidence Interval of \code{n} given \code{alp} along with lower and upper abberation.
#' @return A dataframe with
##'  \item{x }{- Number of successes (positive samples)}
##'  \item{LLR }{ - Likelihood Ratio Lower limit}
##'  \item{ULR }{ - Likelihood Ratio Upper Limit}
##'  \item{LABB }{ - Likelihood Ratio Lower Abberation}
##'  \item{UABB }{ - Likelihood Ratio Upper Abberation}
##'  \item{ZWI }{ - Zero Width Interval}
#' @family Basic methods of CI estimation
#' @seealso \code{\link{prop.test} and \link{binom.test}} for equivalent base Stats R functionality,
#'    \code{\link[binom]{binom.confint}}  provides similar functionality for 11 methods,
#'    \code{\link[PropCIs]{wald2ci}} which provides multiple functions for CI calculation ,
#'    \code{\link[BlakerCI]{binom.blaker.limits}} which calculates Blaker CI which is not covered here and
#'    \code{\link[prevalence]{propCI}} which provides similar functionality.
#' @examples
#' n=5; alp=0.05
#' ciLR(n,alp)
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
ciLR<-function(n,alp)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")

  ####INPUT n
  x=0:n
  k=n+1
  ####INITIALIZATIONS
  mle=0
  cutoff=0
  LLR=0
  ULR=0
  LABB=0
  UABB=0
  ZWI=0

  ###CRITICAL VALUES
  cv=qnorm(1-(alp/2), mean = 0, sd = 1)
  #LIKELIHOOD-RATIO METHOD
  for(i in 1:k)
  {
    likelhd = function(p) dbinom(x[i],n,p)
    loglik = function(p) dbinom(x[i],n,p,log=TRUE)
    mle[i]=optimize(likelhd,c(0,1),maximum=TRUE)$maximum
    cutoff[i]=loglik(mle[i])-(cv^2/2)
    loglik.optim=function(p){abs(cutoff[i]-loglik(p))}
    LLR[i]=optimize(loglik.optim, c(0,mle[i]))$minimum
    ULR[i]=optimize(loglik.optim, c(mle[i],1))$minimum

    if(LLR[i]<0) LABB[i]="YES" else LABB[i]="NO"

    if(ULR[i]>1) UABB[i]="YES" else UABB[i]="NO"

    if(ULR[i]-LLR[i]==0)ZWI[i]="YES" else ZWI[i]="NO"

  }
  return(data.frame(x,LLR,ULR,LABB,UABB,ZWI))
}

##############################################################################
#'  Exact method of CI estimation
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param e - Exact method indicator  in [0, 1] {1: Clopper Pearson, 0.5: Mid P},
#' The input can also be a range of values between 0 and 1.
#' @details  Confidence interval for \code{p} (for all \code{x} = 0, 1, 2 ..\code{n}),
#' based on inverting
#' equal-tailed binomial tests with null hypothesis \deqn{H0: p = p0} and calculated from the
#' cumulative binomial distribution. Exact two sided P-value is usually calculated as
#' \deqn{P= 2[e*Pr(X = x) + min{(Pr(X < x), Pr(X > x))}]} where
#' probabilities are found at null value of p and \eqn{0 <= e <= 1}.
#' The Confidence Interval of \code{n} given \code{alp} along with lower and upper abberation.
#' @return A dataframe with
##'  \item{x }{- Number of successes (positive samples)}
##'  \item{LEX }{ - Exact Lower limit}
##'  \item{UEX }{ - Exact Upper Limit}
##'  \item{LABB }{ - Likelihood Ratio Lower Abberation}
##'  \item{UABB }{ - Likelihood Ratio Upper Abberation}
##'  \item{ZWI }{ - Zero Width Interval}
##'  \item{e }{- Exact method input}
#' @family Basic methods of CI estimation
#' @seealso \code{\link{prop.test} and \link{binom.test}} for equivalent base Stats R functionality,
#'    \code{\link[binom]{binom.confint}}  provides similar functionality for 11 methods,
#'    \code{\link[PropCIs]{wald2ci}} which provides multiple functions for CI calculation ,
#'    \code{\link[BlakerCI]{binom.blaker.limits}} which calculates Blaker CI which is not covered here and
#'    \code{\link[prevalence]{propCI}} which provides similar functionality.
#' @examples
#' n=5; alp=0.05;e=0.5
#' ciEX(n,alp,e) #Mid-p
#' n=5; alp=0.05;e=1 #Clopper-Pearson
#' ciEX(n,alp,e)
#' n=5; alp=0.05;e=c(0.1,0.5,0.95,1) #Range including Mid-p and Clopper-Pearson
#' ciEX(n,alp,e)
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
ciEX=function(n,alp,e)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(e)) stop("'e' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if ((class(e) != "integer") & (class(e) != "numeric") || any(e>1) || any(e<0)) stop("'e' has to be between 0 and 1")
  if (length(e)>10) stop("'e' can have only 10 intervals")

  nvar=length(e)

  res <- data.frame()

  for(i in 1:nvar)
  {
    lu=lufn101(n,alp,e[i])
    res <- rbind(res,lu)
  }
  return(res)
}
lufn101=function(n,alp,e)
{
  x=0:n
  k=n+1
  LEX=0
  UEX=0
  LABB=0
  UABB=0
  ZWI=0

    for(i in 1:k)
    {
    LEX[i]=exlim102l(x[i],n,alp,e)
    UEX[i]=exlim102u(x[i],n,alp,e)
    if(LEX[i]<0) LABB[i]="YES" else LABB[i]="NO"
    if(LEX[i]<0) LEX[i]=0

    if(UEX[i]>1) UABB[i]="YES" else UABB[i]="NO"
    if(UEX[i]>1) UEX[i]=1

    if(UEX[i]-LEX[i]==0)ZWI[i]="YES" else ZWI[i]="NO"
    }
  res=data.frame(x,LEX,UEX,LABB, UABB, ZWI, e)
  return(res)
}

exlim102l=function(x,n,alp,e)
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
exlim102u=function(x,n,alp,e)
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
#'  Wald-T method of CI estimation
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @details  An approximate method based on a t_approximation of the standardized point estimator for all x = 0, 1, 2 ..n; that is the point estimator divided by its estimated standard error.
#' Essential boundary modification is when x = 0 or n, \deqn{\hat{p}=\frac{(x+2)}{(n+4)}}
#' @return A dataframe with
##'  \item{x }{- Number of successes (positive samples)}
##'  \item{LTW }{ - Wald-T Lower limit}
##'  \item{UTW }{ - Wald-T Upper Limit}
##'  \item{LABB }{ - Wald-T Lower Abberation}
##'  \item{UABB }{ - Wald-T Upper Abberation}
##'  \item{ZWI }{ - Zero Width Interval}
#' @family Basic methods of CI estimation
#' @seealso \code{\link{prop.test} and \link{binom.test}} for equivalent base Stats R functionality,
#'    \code{\link[binom]{binom.confint}}  provides similar functionality for 11 methods,
#'    \code{\link[PropCIs]{wald2ci}} which provides multiple functions for CI calculation ,
#'    \code{\link[BlakerCI]{binom.blaker.limits}} which calculates Blaker CI which is not covered here and
#'    \code{\link[prevalence]{propCI}} which provides similar functionality.
#' @examples
#' n=5; alp=0.05
#' ciTW(n,alp)
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
ciTW<-function(n,alp)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")

  ####INPUT n
  x=0:n
  k=n+1
  ####INITIALIZATIONS
  pTW=0
  qTW=0
  seTW=0
  LTW=0
  UTW=0
  DOF=0
  cv=0
  LABB=0
  UABB=0
  ZWI=0

  #MODIFIED_t-WALD METHOD
  for(i in 1:k)
  {
    if(x[i]==0||x[i]==n)
    {
      pTW[i]=(x[i]+2)/(n+4)
      qTW[i]=1-pTW[i]
    }else
    {
      pTW[i]=x[i]/n
      qTW[i]=1-pTW[i]
    }
    f1=function(p,n) p*(1-p)/n
    f2=function(p,n) (p*(1-p)/(n^3))+(p+((6*n)-7)*(p^2)+(4*(n-1)*(n-3)*(p^3))-(2*(n-1)*((2*n)-3)*(p^4)))/(n^5)-(2*(p+((2*n)-3)*(p^2)-2*(n-1)*(p^3)))/(n^4)
    DOF[i]=2*((f1(pTW[i],n))^2)/f2(pTW[i],n)
    cv[i]=qt(1-(alp/2), df=DOF[i])
    seTW[i]=cv[i]*sqrt(f1(pTW[i],n))
    LTW[i]=pTW[i]-(seTW[i])
    UTW[i]=pTW[i]+(seTW[i])

    if(LTW[i]<0) LABB[i]="YES" else LABB[i]="NO"
    if(LTW[i]<0) LTW[i]=0

    if(UTW[i]>1) UABB[i]="YES" else  UABB[i]="NO"
    if(UTW[i]>1) UTW[i]=1

    if(UTW[i]-LTW[i]==0)ZWI[i]="YES" else ZWI[i]="NO"
  }
  return(data.frame(x,LTW,UTW,LABB,UABB,ZWI))
}
#####################################################################################


###############
#'  Logit Wald method of CI estimation
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @details  The Confidence Interval of \code{n} given \code{alp} along with lower and upper abberation.
#' @return A dataframe with
##'  \item{x }{- Number of successes (positive samples)}
##'  \item{LLT }{ - Logit Wald Lower limit}
##'  \item{ULT }{ - Logit Wald Upper Limit}
##'  \item{LABB }{ - Logit Wald Lower Abberation}
##'  \item{UABB }{ - Logit Wald Upper Abberation}
##'  \item{ZWI }{ - Zero Width Interval}
#' @family Basic methods of CI estimation
#' @seealso \code{\link{prop.test} and \link{binom.test}} for equivalent base Stats R functionality,
#'    \code{\link[binom]{binom.confint}}  provides similar functionality for 11 methods,
#'    \code{\link[PropCIs]{wald2ci}} which provides multiple functions for CI calculation ,
#'    \code{\link[BlakerCI]{binom.blaker.limits}} which calculates Blaker CI which is not covered here and
#'    \code{\link[prevalence]{propCI}} which provides similar functionality.
#' @examples
#' n=5; alp=0.05
#' ciLT(n,alp)
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
ciLT<-function(n,alp)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")

  ####INPUT n
  x=0:n
  k=n+1
  ####INITIALIZATIONS
  pLT=0
  qLT=0
  seLT=0
  lgit=0
  LLT=0
  ULT=0
  LABB=0
  UABB=0
  ZWI=0
  ###CRITICAL VALUES
  cv=qnorm(1-(alp/2), mean = 0, sd = 1)
  #LOGIT-WALD METHOD
  pLT[1]=0
  qLT[1]=1
  LLT[1] = 0
  ULT[1] = 1-((alp/2)^(1/n))

  pLT[k]=1
  qLT[k]=0
  LLT[k]= (alp/2)^(1/n)
  ULT[k]=1

  for(j in 1:(k-2))
  {
    pLT[j+1]=x[j+1]/n
    qLT[j+1]=1-pLT[j+1]
    lgit[j+1]=log(pLT[j+1]/qLT[j+1])
    seLT[j+1]=sqrt(pLT[j+1]*qLT[j+1]*n)
    LLT[j+1]=1/(1+exp(-lgit[j+1]+(cv/seLT[j+1])))
    ULT[j+1]=1/(1+exp(-lgit[j+1]-(cv/seLT[j+1])))
  }
  for(i in 1:k)
  {
    if(LLT[i]<0) LABB[i]="YES" else LABB[i]="NO"
    if(LLT[i]<0) LLT[i]=0

    if(ULT[i]>1) UABB[i]="YES" else UABB[i]="NO"
    if(ULT[i]>1) ULT[i]=1

    if(ULT[i]-LLT[i]==0)ZWI[i]="YES" else ZWI[i]="NO"
  }
  return(data.frame(x,LLT,ULT,LABB,UABB,ZWI))
}
#####################################################################
#' Bayesian method of CI estimation with different or same parameteric values
#' for Beta prior distribution
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Shape parameter 1 for prior Beta distribution in Bayesian model.
#' Can also be a vector of length n+1 priors.
#' @param b - Shape parameter 2 for prior Beta distribution in Bayesian model.
#' Can also be a vector of length n+1 priors.
#' @details  Highest Probability Density (HPD) and two tailed intervals are provided for all
#' \eqn{xi = 0, 1, 2 ..n} based on the conjugate prior \eqn{\beta(ai, bi) (i = 1, 2..n+1)}
#'  for the probability of success \code{p} of the binomial distribution so that the posterior
#'  is \eqn{\beta(xi + ai, n - xi + bi)}.
#' @return A dataframe with
##'  \item{x }{- Number of successes (positive samples)}
##'  \item{pomean }{ - Posterior mean}
##'  \item{LBAQ }{ - Lower limits of Quantile based intervals}
##'  \item{UBAQ }{ - Upper limits of Quantile based intervals}
##'  \item{LBAH }{ - Lower limits of HPD intervals}
##'  \item{UBAH }{ - Upper limits of HPD intervals}
#' @family Basic methods of CI estimation
#' @seealso \code{\link{prop.test} and \link{binom.test}} for equivalent base Stats R functionality,
#'    \code{\link[binom]{binom.confint}}  provides similar functionality for 11 methods,
#'    \code{\link[PropCIs]{wald2ci}} which provides multiple functions for CI calculation ,
#'    \code{\link[BlakerCI]{binom.blaker.limits}} which calculates Blaker CI which is not covered here and
#'    \code{\link[prevalence]{propCI}} which provides similar functionality.
#' @examples
#' n=5; alp=0.05; a=0.5;b=0.5;
#' ciBA(n,alp,a,b)
#' n=5; alp=0.05; a=c(0.5,2,1,1,2,0.5);b=c(0.5,2,1,1,2,0.5)
#' ciBA(n,alp,a,b)
#' @references
#' [1] 2002 Gelman A, Carlin  JB, Stern HS and Dunson DB
#' Bayesian Data Analysis, Chapman & Hall/CRC
#' [2] 2006  Ghosh M, Delampady M and Samanta T.
#' An introduction to Bayesian analysis: Theory and Methods. Springer, New York
#' @export
ciBA<-function(n,alp,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")

  if(length(a)>1 || length(b)>1){
  BAdf=ciBAD(n,alp,a,b)}
 else{
  if ((class(a) != "integer") & (class(a) != "numeric") || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || b<0  ) stop("'b' has to be greater than or equal to 0")

  ####INPUT n
  x=0:n
  k=n+1
  ####INITIALIZATIONS
  LBAQ=0
  UBAQ=0
  LBAH=0
  UBAH=0
  pomean=0
  ##############
  #  library(TeachingDemos)				#To get HPDs
  for(i in 1:k)
  {
    if (!requireNamespace("TeachingDemos", quietly = TRUE)) {
      stop("TeachingDemos needed for this function to work. Please install it.",
           call. = FALSE)}
    #Posterior mean
    pomean[i]=(x[i]+a)/(n+a+b)

    #Quantile Based Intervals
    LBAQ[i]=qbeta(alp/2,x[i]+a,n-x[i]+b)
    UBAQ[i]=qbeta(1-(alp/2),x[i]+a,n-x[i]+b)

    LBAH[i]=TeachingDemos::hpd(qbeta,shape1=x[i]+a,shape2=n-x[i]+b,conf=1-alp)[1]
    UBAH[i]=TeachingDemos::hpd(qbeta,shape1=x[i]+a,shape2=n-x[i]+b,conf=1-alp)[2]

  }
    BAdf=data.frame(x,pomean,LBAQ,UBAQ,LBAH,UBAH)
 }
    return(BAdf)
}

#######################################################################
ciBAD<-function(n,alp,a,b)
{
  if ((class(a) != "integer") & (class(a) != "numeric") || any(a < 0)) stop("'a' has to be a set of positive numeric vectors")
  if ((class(b) != "integer") & (class(b) != "numeric") || any(b < 0)) stop("'b' has to be a set of positive numeric vectors")
  if (length(a) <  n || length(b) <  n ) stop("'a' and 'b' vectors have to be equal to length n")


  ####INPUT n
  x=0:n
  k=n+1
  ####INITIALIZATIONS
  LBAQ=0
  UBAQ=0
  LBAH=0
  UBAH=0
  pomean=0
  ##############
  #  library(TeachingDemos)				#To get HPDs
  for(i in 1:k)
  {
    #Posterior mean
    pomean[i]=(x[i]+a[i])/(n+a[i]+b[i])

    #Quantile Based Intervals
    LBAQ[i]=qbeta(alp/2,x[i]+a[i],n-x[i]+b[i])
    UBAQ[i]=qbeta(1-(alp/2),x[i]+a[i],n-x[i]+b[i])

    LBAH[i]=TeachingDemos::hpd(qbeta,shape1=x[i]+a[i],shape2=n-x[i]+b[i],conf=1-alp)[1]
    UBAH[i]=TeachingDemos::hpd(qbeta,shape1=x[i]+a[i],shape2=n-x[i]+b[i],conf=1-alp)[2]

  }
  return(data.frame(x,pomean,LBAQ,UBAQ,LBAH,UBAH))
}

#####################################################################################
#' CI estimation of 6 base methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @details  The Confidence Interval of 6 base methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine) for \code{n} given \code{alp}.
#' @return A dataframe with
##'  \item{method }{- Name of the method}
##'  \item{x }{- Number of successes (positive samples)}
##'  \item{LLT }{ - Lower limit}
##'  \item{ULT }{ - Upper Limit}
##'  \item{LABB }{ - Lower Abberation}
##'  \item{UABB }{ - Upper Abberation}
##'  \item{ZWI }{ - Zero Width Interval}
#' @family Basic methods of CI estimation
#' @seealso \code{\link{prop.test} and \link{binom.test}} for equivalent base Stats R functionality,
#'    \code{\link[binom]{binom.confint}}  provides similar functionality for 11 methods,
#'    \code{\link[PropCIs]{wald2ci}} which provides multiple functions for CI calculation ,
#'    \code{\link[BlakerCI]{binom.blaker.limits}} which calculates Blaker CI which is not covered here and
#'    \code{\link[prevalence]{propCI}} which provides similar functionality.
#' @examples
#' n=5; alp=0.05;
#' ciAll(n,alp)
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
ciAll<-function(n,alp)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")

  #### Calling functions and creating df
  WaldCI.df = ciWD(n,alp)
  ArcSineCI.df = ciAS(n,alp)
  LRCI.df = ciLR(n,alp)
  ScoreCI.df = ciSC(n,alp)
  WaldTCI.df = ciTW(n,alp)
  LogitWald.df = ciLT(n,alp)

  WaldCI.df$method = as.factor("Wald")
  ArcSineCI.df$method = as.factor("ArcSine")
  LRCI.df$method = as.factor("Likelihood")
  ScoreCI.df$method = as.factor("Score")
  WaldTCI.df$method = as.factor("Wald-T")
  LogitWald.df$method = as.factor("Logit-Wald")

  Generic.1 = data.frame(method = WaldCI.df$method, x=WaldCI.df$x, LowerLimit = WaldCI.df$LWD, UpperLimit = WaldCI.df$UWD, LowerAbb = WaldCI.df$LABB, UpperAbb = WaldCI.df$UABB, ZWI = WaldCI.df$ZWI)
  Generic.2 = data.frame(method = ArcSineCI.df$method, x=ArcSineCI.df$x, LowerLimit = ArcSineCI.df$LAS, UpperLimit = ArcSineCI.df$UAS, LowerAbb = ArcSineCI.df$LABB, UpperAbb = ArcSineCI.df$UABB, ZWI = ArcSineCI.df$ZWI)
  Generic.3 = data.frame(method = LRCI.df$method, x=LRCI.df$x, LowerLimit = LRCI.df$LLR, UpperLimit = LRCI.df$ULR, LowerAbb = LRCI.df$LABB, UpperAbb = LRCI.df$UABB, ZWI = LRCI.df$ZWI)
  Generic.4 = data.frame(method = ScoreCI.df$method, x=ScoreCI.df$x, LowerLimit = ScoreCI.df$LSC, UpperLimit = ScoreCI.df$USC, LowerAbb = ScoreCI.df$LABB, UpperAbb = ScoreCI.df$UABB, ZWI = ScoreCI.df$ZWI)
  Generic.5 = data.frame(method = WaldTCI.df$method, x=WaldTCI.df$x, LowerLimit = WaldTCI.df$LTW, UpperLimit = WaldTCI.df$UTW, LowerAbb = WaldTCI.df$LABB, UpperAbb = WaldTCI.df$UABB, ZWI = WaldTCI.df$ZWI)
  Generic.6 = data.frame(method = LogitWald.df$method, x=LogitWald.df$x, LowerLimit = LogitWald.df$LLT, UpperLimit = LogitWald.df$ULT, LowerAbb = LogitWald.df$LABB, UpperAbb = LogitWald.df$UABB, ZWI = LogitWald.df$ZWI)

  Final.df= rbind(Generic.1,Generic.2,Generic.3,Generic.4,Generic.5, Generic.6)

  return(Final.df)
}
