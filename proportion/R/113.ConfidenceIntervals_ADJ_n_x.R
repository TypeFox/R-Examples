#' Adjusted Wald method of CI estimation
#' @param x - Number of successes
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @details  Given data \code{x} and \code{n} are modified as \eqn{x + h} and \eqn{n + (2*h)}
#'  respectively, where \eqn{h > 0} then Wald-type interval is applied for the given \code{x} and \code{n}.
#' @return A dataframe with
#'  \item{x}{  Number of successes (positive samples)}
#'  \item{LAWDx }{   Adjusted Wald Lower limit}
#'  \item{UAWDx }{   Adjusted Wald Upper Limit}
#'  \item{LABB }{   Adjusted Wald Lower Abberation}
#'  \item{UABB }{   Adjusted Wald Upper Abberation}
#'  \item{ZWI }{   Zero Width Interval}
#' @family Adjusted methods of CI estimation  given x & n
#' @seealso \code{\link{prop.test} and \link{binom.test}} for equivalent base Stats R functionality,
#'    \code{\link[binom]{binom.confint}}  provides similar functionality for 11 methods,
#'    \code{\link[PropCIs]{wald2ci}} which provides multiple functions for CI calculation ,
#'    \code{\link[BlakerCI]{binom.blaker.limits}} which calculates Blaker CI which is not covered here and
#'    \code{\link[prevalence]{propCI}} which provides similar functionality.
#' @examples
#' x= 5; n=5; alp=0.05; h=2
#' ciAWDx(x,n,alp,h)
#' @references
#' [1] 1998 Agresti A and Coull BA.
#' Approximate is better than "Exact" for interval estimation of binomial proportions.
#' The American Statistician: 52; 119 - 126.
#'
#' [2] 1998 Newcombe RG.
#' Two-sided confidence intervals for the single proportion: Comparison of seven methods.
#' Statistics in Medicine: 17; 857 - 872.
#'
#' [3] 2008 Pires, A.M., Amado, C.
#' Interval Estimators for a Binomial Proportion: Comparison of Twenty Methods.
#' REVSTAT - Statistical Journal, 6, 165-197.
#' @export
#1.WALD
ciAWDx<-function(x,n,alp,h)
{
  if (missing(x)) stop("'x' is missing")
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (((class(x) != "integer") & (class(x) != "numeric")) || (x<0) || x>n || length(x)>1) stop("'x' has to be a positive integer between 0 and n")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h) >1|| h<0  ) stop("'h' has to be greater than or equal to 0")

###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#WALD METHOD
y=x+h
m=n+(2*h)
pAW=y/m
qAW=1-pAW
seAW=sqrt(pAW*qAW/m)
LAWDx=pAW-(cv*seAW)
UAWDx=pAW+(cv*seAW)
if(LAWDx<0) LABB="YES" else LABB="NO"
if(LAWDx<0) LAWDx=0
if(UAWDx>1) UABB="YES" else  UABB="NO"
if(UAWDx>1) UAWDx=1
if(UAWDx-LAWDx==0)ZWI="YES" else ZWI="NO"
return(data.frame(x,LAWDx,UAWDx,LABB,UABB,ZWI))
}
##############################################################################
#' Adjusted Score method of CI estimation
#' @param x - Number of successes
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @details  A score test approach is used after the given data \code{x} and \code{n} are
#' modified as \eqn{x + h} and \eqn{n + (2*h)} respectively, where \eqn{h > 0} and for the
#'  given \code{x} and \code{n}.
#' @return A dataframe with
#'  \item{x}{  Number of successes (positive samples)}
#'  \item{LASCx }{   Score Lower limit}
#'  \item{UASCx }{   Score Upper Limit}
#'  \item{LABB }{   Score Lower Abberation}
#'  \item{UABB }{   Score Upper Abberation}
#'  \item{ZWI }{   Zero Width Interval}
#' @family Adjusted methods of CI estimation  given x & n
#' @seealso \code{\link{prop.test} and \link{binom.test}} for equivalent base Stats R functionality,
#'    \code{\link[binom]{binom.confint}}  provides similar functionality for 11 methods,
#'    \code{\link[PropCIs]{wald2ci}} which provides multiple functions for CI calculation ,
#'    \code{\link[BlakerCI]{binom.blaker.limits}} which calculates Blaker CI which is not covered here and
#'    \code{\link[prevalence]{propCI}} which provides similar functionality.
#' @examples
#' x=5; n=5; alp=0.05; h=2
#' ciASCx(x,n,alp,h)
#' @references
#' [1] 1998 Agresti A and Coull BA.
#' Approximate is better than "Exact" for interval estimation of binomial proportions.
#' The American Statistician: 52; 119 - 126.
#'
#' [2] 1998 Newcombe RG.
#' Two-sided confidence intervals for the single proportion: Comparison of seven methods.
#' Statistics in Medicine: 17; 857 - 872.
#'
#' [3] 2008 Pires, A.M., Amado, C.
#' Interval Estimators for a Binomial Proportion: Comparison of Twenty Methods.
#' REVSTAT - Statistical Journal, 6, 165-197.
#' @export
#2.SCORE
ciASCx<-function(x,n,alp,h)
{
  if (missing(x)) stop("'x' is missing")
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (((class(x) != "integer") & (class(x) != "numeric")) || (x<0) || x>n || length(x)>1) stop("'x' has to be a positive integer between 0 and n")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h) >1|| h<0  ) stop("'h' has to be greater than or equal to 0")

y=x+h
m=n+(2*h)

###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
cv1=(cv^2)/(2*m)
cv2=(cv/(2*m))^2

#AASCORE (WILSON) METHOD
pAS=y/m
qAS=1-pAS
seAS=sqrt((pAS*qAS/m)+cv2)
LASCx=(m/(m+(cv)^2))*((pAS+cv1)-(cv*seAS))
UASCx=(m/(m+(cv)^2))*((pAS+cv1)+(cv*seAS))
if(LASCx<0) LABB="YES" else LABB="NO"
if(LASCx<0) LASCx=0

if(UASCx>1) UABB="YES" else  UABB="NO"
if(UASCx>1) UASCx=1

if(UASCx-LASCx==0)ZWI="YES" else ZWI="NO"

return(data.frame(x,LASCx,UASCx,LABB,UABB,ZWI))
}
##############################################################################
#' Adjusted ArcSine method of CI estimation
#' @param x - Number of successes
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @details  Wald-type interval for the arcsine transformation of the parameter \code{p}
#' for the modified data \eqn{x + h} and \eqn{n + (2*h)} , where \eqn{h > 0} and for
#' the given \code{x} and \code{n}.
#' @return A dataframe with
#'  \item{x}{  Number of successes (positive samples)}
#'  \item{LAASx }{   ArcSine Lower limit}
#'  \item{UAASx }{   ArcSine Upper Limit}
#'  \item{LABB }{   ArcSine Lower Abberation}
#'  \item{UABB }{   ArcSine Upper Abberation}
#'  \item{ZWI }{   Zero Width Interval}
#' @family Adjusted methods of CI estimation  given x & n
#' @seealso \code{\link{prop.test} and \link{binom.test}} for equivalent base Stats R functionality,
#'    \code{\link[binom]{binom.confint}}  provides similar functionality for 11 methods,
#'    \code{\link[PropCIs]{wald2ci}} which provides multiple functions for CI calculation ,
#'    \code{\link[BlakerCI]{binom.blaker.limits}} which calculates Blaker CI which is not covered here and
#'    \code{\link[prevalence]{propCI}} which provides similar functionality.
#' @examples
#' x=5; n=5; alp=0.05;h=2
#' ciAASx(x,n,alp,h)
#' @references
#' [1] 1998 Agresti A and Coull BA.
#' Approximate is better than "Exact" for interval estimation of binomial proportions.
#' The American Statistician: 52; 119 - 126.
#'
#' [2] 1998 Newcombe RG.
#' Two-sided confidence intervals for the single proportion: Comparison of seven methods.
#' Statistics in Medicine: 17; 857 - 872.
#'
#' [3] 2008 Pires, A.M., Amado, C.
#' Interval Estimators for a Binomial Proportion: Comparison of Twenty Methods.
#' REVSTAT - Statistical Journal, 6, 165-197.
#' @export
#3.ARCSINE
ciAASx<-function(x,n,alp,h)
{
  if (missing(x)) stop("'x' is missing")
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (((class(x) != "integer") & (class(x) != "numeric")) || (x<0) || x>n || length(x)>1) stop("'x' has to be a positive integer between 0 and n")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h) >1|| h<0  ) stop("'h' has to be greater than or equal to 0")

y=x+h
m=n+(2*h)
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#ARC-SINE METHOD
pA=y/m
#qA=1-pA
seA=cv/sqrt(4*m)
LAASx=(sin(asin(sqrt(pA))-seA))^2
UAASx=(sin(asin(sqrt(pA))+seA))^2

if(LAASx<0) LABB="YES" else LABB="NO"
if(LAASx<0) LAASx=0

if(UAASx>1) UABB="YES" else  UABB="NO"
if(UAASx>1) UAASx=1

if(UAASx-LAASx==0)ZWI="YES" else ZWI="NO"
return(data.frame(x,LAASx,UAASx,LABB,UABB,ZWI))
}
##############################################################################
#' AdjustedLikelyhood Ratio method of CI estimation
#' @param x - Number of successes
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @details  Likelihood ratio limits for the data \eqn{x + h} and \eqn{n + (2*h)}
#' instead of the given \code{x} and \code{n}, where \code{h} is a positive integer
#' \eqn{(1, 2.)} and for the given \code{x} and \code{n}.
#' @return A dataframe with
#'  \item{x}{  Number of successes (positive samples)}
#'  \item{LALRx }{   Likelyhood Ratio Lower limit}
#'  \item{UALRx }{   Likelyhood Ratio Upper Limit}
#'  \item{LABB }{   Likelyhood Ratio Lower Abberation}
#'  \item{UABB }{   Likelyhood Ratio Upper Abberation}
#'  \item{ZWI }{   Zero Width Interval}
#' @family Adjusted methods of CI estimation  given x & n
#' @seealso \code{\link{prop.test} and \link{binom.test}} for equivalent base Stats R functionality,
#'    \code{\link[binom]{binom.confint}}  provides similar functionality for 11 methods,
#'    \code{\link[PropCIs]{wald2ci}} which provides multiple functions for CI calculation ,
#'    \code{\link[BlakerCI]{binom.blaker.limits}} which calculates Blaker CI which is not covered here and
#'    \code{\link[prevalence]{propCI}} which provides similar functionality.
#' @examples
#' x=5; n=5; alp=0.05;h=2
#' ciALRx(x,n,alp,h)
#' @references
#' [1] 1998 Agresti A and Coull BA.
#' Approximate is better than "Exact" for interval estimation of binomial proportions.
#' The American Statistician: 52; 119 - 126.
#'
#' [2] 1998 Newcombe RG.
#' Two-sided confidence intervals for the single proportion: Comparison of seven methods.
#' Statistics in Medicine: 17; 857 - 872.
#'
#' [3] 2008 Pires, A.M., Amado, C.
#' Interval Estimators for a Binomial Proportion: Comparison of Twenty Methods.
#' REVSTAT - Statistical Journal, 6, 165-197.
#' @export
#4.LIKELIHOOD RATIO
ciALRx<-function(x,n,alp,h)
{
  if (missing(x)) stop("'x' is missing")
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (((class(x) != "integer") & (class(x) != "numeric")) || (x<0) || x>n || length(x)>1) stop("'x' has to be a positive integer between 0 and n")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h) >1|| h<0  || !(h%%1 ==0)) stop("'h' has to be an integer greater than or equal to 0")

####INPUT
y=x
y1=y+h
n1=n+(2*h)

###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#LIKELIHOOD-RATIO METHOD
likelhd = function(p) dbinom(y1,n1,p)
loglik = function(p) dbinom(y1,n1,p,log=TRUE)
mle=optimize(likelhd,c(0,1),maximum=TRUE)$maximum
cutoff=loglik(mle)-(cv^2/2)
loglik.optim=function(p){abs(cutoff-loglik(p))}
LALRx=optimize(loglik.optim, c(0,mle))$minimum
UALRx=optimize(loglik.optim, c(mle,1))$minimum

if(LALRx<0) LABB="YES" else LABB="NO"
if(LALRx<0) LALRx=0

if(UALRx>1) UABB="YES" else  UABB="NO"
if(UALRx>1) UALRx=1

if(UALRx-LALRx==0)ZWI="YES" else ZWI="NO"
return(data.frame(x,LALRx,UALRx,LABB,UABB,ZWI))
}
######################################################################
#' AdjustedWALD-T method of CI estimation
#' @param x - Number of successes
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @details  Given data \code{x} and \code{n} are modified as \eqn{x + h} and \eqn{n + (2*h)}
#'  respectively, where \eqn{h > 0} then approximate method based on a t_approximation of the
#'  standardized point estimator for the given \code{x} and \code{n}.
#' @return A dataframe with
#'  \item{x}{  Number of successes (positive samples)}
#'  \item{LATWx }{   T-Wald Lower limit}
#'  \item{UATWx }{   T-Wald Upper Limit}
#'  \item{LABB }{   T-Wald Lower Abberation}
#'  \item{UABB }{   T-Wald Upper Abberation}
#'  \item{ZWI }{   Zero Width Interval}
#' @family Adjusted methods of CI estimation  given x & n
#' @seealso \code{\link{prop.test} and \link{binom.test}} for equivalent base Stats R functionality,
#'    \code{\link[binom]{binom.confint}}  provides similar functionality for 11 methods,
#'    \code{\link[PropCIs]{wald2ci}} which provides multiple functions for CI calculation ,
#'    \code{\link[BlakerCI]{binom.blaker.limits}} which calculates Blaker CI which is not covered here and
#'    \code{\link[prevalence]{propCI}} which provides similar functionality.
#' @examples
#' x=5; n=5; alp=0.05;h=2
#' ciATWx(x,n,alp,h)
#' @references
#' [1] 1998 Agresti A and Coull BA.
#' Approximate is better than "Exact" for interval estimation of binomial proportions.
#' The American Statistician: 52; 119 - 126.
#'
#' [2] 1998 Newcombe RG.
#' Two-sided confidence intervals for the single proportion: Comparison of seven methods.
#' Statistics in Medicine: 17; 857 - 872.
#'
#' [3] 2008 Pires, A.M., Amado, C.
#' Interval Estimators for a Binomial Proportion: Comparison of Twenty Methods.
#' REVSTAT - Statistical Journal, 6, 165-197.
#' @export
#5.ADJUSTED WALD-T
ciATWx<-function(x,n,alp,h)
{
  if (missing(x)) stop("'x' is missing")
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (((class(x) != "integer") & (class(x) != "numeric")) || (x<0) || x>n || length(x)>1) stop("'x' has to be a positive integer between 0 and n")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h) >1|| h<0  ) stop("'h' has to be greater than or equal to 0")

y=x+h
n1=n+(2*h)
								#Coverage probabilty
#MODIFIED_t-ADJ_WALD METHOD
pATWx=y/n1
#qATWx=1-pATWx
f1=function(p,n) p*(1-p)/n
f2=function(p,n) (p*(1-p)/(n^3))+(p+((6*n)-7)*(p^2)+(4*(n-1)*(n-3)*(p^3))-(2*(n-1)*((2*n)-3)*(p^4)))/(n^5)-(2*(p+((2*n)-3)*(p^2)-2*(n-1)*(p^3)))/(n^4)
DOFx=2*((f1(pATWx,n1))^2)/f2(pATWx,n1)
cvx=qt(1-(alp/2), df=DOFx)
seATWx=cvx*sqrt(f1(pATWx,n1))
LATWx=pATWx-(seATWx)
UATWx=pATWx+(seATWx)

if(LATWx<0) LABB="YES" else LABB="NO"
if(LATWx<0) LATWx=0

if(UATWx>1) UABB="YES" else  UABB="NO"
if(UATWx>1) UATWx=1

if(UATWx-LATWx==0)ZWI="YES" else ZWI="NO"

return(data.frame(x,LATWx,UATWx,LABB,UABB,ZWI))
}
#####################################################################
#' Adjusted Logit-Wald method of CI estimation
#' @param x - Number of successes
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @details  Wald-type interval for the logit transformation \eqn{log(p/1-p)}
#' of the parameter p for the modified data eqn{x + h} and \eqn{n + (2*h)} ,
#' where \eqn{h > 0}  and the given \code{x} and \code{n}.
#' @return A dataframe with
#'  \item{x}{  Number of successes (positive samples)}
#'  \item{LALTx }{   Logit Wald Lower limit}
#'  \item{UALTx }{   Logit Wald Upper Limit}
#'  \item{LABB }{   Logit Wald Lower Abberation}
#'  \item{UABB }{   Logit Wald Upper Abberation}
#'  \item{ZWI }{   Zero Width Interval}
#' @family Adjusted methods of CI estimation  given x & n
#' @seealso \code{\link{prop.test} and \link{binom.test}} for equivalent base Stats R functionality,
#'    \code{\link[binom]{binom.confint}}  provides similar functionality for 11 methods,
#'    \code{\link[PropCIs]{wald2ci}} which provides multiple functions for CI calculation ,
#'    \code{\link[BlakerCI]{binom.blaker.limits}} which calculates Blaker CI which is not covered here and
#'    \code{\link[prevalence]{propCI}} which provides similar functionality.
#' @examples
#' x=5; n=5; alp=0.05;h=2
#' ciALTx(x,n,alp,h)
#' @references
#' [1] 1998 Agresti A and Coull BA.
#' Approximate is better than "Exact" for interval estimation of binomial proportions.
#' The American Statistician: 52; 119 - 126.
#'
#' [2] 1998 Newcombe RG.
#' Two-sided confidence intervals for the single proportion: Comparison of seven methods.
#' Statistics in Medicine: 17; 857 - 872.
#'
#' [3] 2008 Pires, A.M., Amado, C.
#' Interval Estimators for a Binomial Proportion: Comparison of Twenty Methods.
#' REVSTAT - Statistical Journal, 6, 165-197.
#' @export
#6.ADJUSTED LOGIT-WALD
ciALTx<-function(x,n,alp,h)
{
  if (missing(x)) stop("'x' is missing")
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (((class(x) != "integer") & (class(x) != "numeric")) || (x<0) || x>n || length(x)>1) stop("'x' has to be a positive integer between 0 and n")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h) >1|| h<0  ) stop("'h' has to be greater than or equal to 0")

y=x+h
n1=n+(2*h)

###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#LOGIT-WALD METHOD
pALTx=y/n1
qALTx=1-pALTx
lgitx=log(pALTx/qALTx)
seALTx=sqrt(pALTx*qALTx*n1)
LALTx=1/(1+exp(-lgitx+(cv/seALTx)))
UALTx=1/(1+exp(-lgitx-(cv/seALTx)))

if(LALTx<0) LABB="YES" else LABB="NO"
if(LALTx<0) LALTx=0

if(UALTx>1) UABB="YES" else  UABB="NO"
if(UALTx>1) UALTx=1

if(UALTx-LALTx==0)ZWI="YES" else ZWI="NO"

return(data.frame(x,LALTx,UALTx,LABB,UABB,ZWI))
}

#####################################################################
#' CI estimation of 6 adjusted methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine) given x and n
#' @param x - Number of success
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @details  The Confidence Interval of using 6 adjusted methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine) for \code{n} given \code{alp}, \code{x} and \code{h}
#' @return A dataframe with
#'  \item{name }{- Name of the method}
#'  \item{x }{- Number of successes (positive samples)}
#'  \item{LLT }{ - Lower limit}
#'  \item{ULT }{ - Upper Limit}
#'  \item{LABB }{ - Lower Abberation}
#'  \item{UABB }{ - Upper Abberation}
#'  \item{ZWI }{ - Zero Width Interval}
#' @family Adjusted methods of CI estimation  given x & n
#' @seealso \code{\link{prop.test} and \link{binom.test}} for equivalent base Stats R functionality,
#'    \code{\link[binom]{binom.confint}}  provides similar functionality for 11 methods,
#'    \code{\link[PropCIs]{wald2ci}} which provides multiple functions for CI calculation ,
#'    \code{\link[BlakerCI]{binom.blaker.limits}} which calculates Blaker CI which is not covered here and
#'    \code{\link[prevalence]{propCI}} which provides similar functionality.
#' @examples
#' x=5; n=5; alp=0.05;h=2
#' ciAAllx(x,n,alp,h)
#' @references
#' [1] 1998 Agresti A and Coull BA.
#' Approximate is better than "Exact" for interval estimation of binomial proportions.
#' The American Statistician: 52; 119 - 126.
#'
#' [2] 1998 Newcombe RG.
#' Two-sided confidence intervals for the single proportion: Comparison of seven methods.
#' Statistics in Medicine: 17; 857 - 872.
#'
#' [3] 2008 Pires, A.M., Amado, C.
#' Interval Estimators for a Binomial Proportion: Comparison of Twenty Methods.
#' REVSTAT - Statistical Journal, 6, 165-197.
#' @export
#9.All methods
ciAAllx<-function(x,n,alp,h)
{
  if (missing(x)) stop("'x' is missing")
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (((class(x) != "integer") & (class(x) != "numeric")) || (x<0) || x>n || length(x)>1) stop("'x' has to be a positive integer between 0 and n")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h) >1|| h<0  || !(h%%1 ==0)) stop("'h' has to be an integer greater than or equal to 0")

  #### Calling functions and creating df
  WaldCI.df    = ciAWDx(x,n,alp,h)
  ArcSineCI.df = ciAASx(x,n,alp,h)
  LRCI.df      = ciALRx(x,n,alp,h)
  ScoreCI.df   = ciASCx(x,n,alp,h)
  WaldLCI.df   = ciALTx(x,n,alp,h)
  AdWaldCI.df  = ciATWx(x,n,alp,h)

  WaldCI.df$method    = as.factor("Adj-Wald")
  ArcSineCI.df$method = as.factor("Adj-ArcSine")
  LRCI.df$method      = as.factor("Adj-Liklihood")
  WaldLCI.df$method    = as.factor("Adj-Logit Wald")
  ScoreCI.df$method   = as.factor("Adj-Score")
  AdWaldCI.df$method  = as.factor("Adj-Wald-T")

  Generic.1 = data.frame(method = WaldCI.df$method, x=WaldCI.df$x, LowerLimit = WaldCI.df$LAWDx, UpperLimit = WaldCI.df$UAWDx, LowerAbb = WaldCI.df$LABB, UpperAbb = WaldCI.df$UABB, ZWI = WaldCI.df$ZWI)
  Generic.2 = data.frame(method = ArcSineCI.df$method, x=ArcSineCI.df$x, LowerLimit = ArcSineCI.df$LAASx, UpperLimit = ArcSineCI.df$UAASx, LowerAbb = ArcSineCI.df$LABB, UpperAbb = ArcSineCI.df$UABB, ZWI = ArcSineCI.df$ZWI)
  Generic.3 = data.frame(method = LRCI.df$method, x=LRCI.df$x, LowerLimit = LRCI.df$LALRx, UpperLimit = LRCI.df$UALRx, LowerAbb = LRCI.df$LABB, UpperAbb = LRCI.df$UABB, ZWI = LRCI.df$ZWI)
  Generic.4 = data.frame(method = ScoreCI.df$method, x=ScoreCI.df$x, LowerLimit = ScoreCI.df$LASCx, UpperLimit = ScoreCI.df$UASCx, LowerAbb = ScoreCI.df$LABB, UpperAbb = ScoreCI.df$UABB, ZWI = ScoreCI.df$ZWI)
  Generic.5 = data.frame(method = WaldLCI.df$method, x=WaldLCI.df$x, LowerLimit = WaldLCI.df$LALTx, UpperLimit = WaldLCI.df$UALTx, LowerAbb = WaldLCI.df$LABB, UpperAbb = WaldLCI.df$UABB, ZWI = WaldLCI.df$ZWI)
  Generic.6 = data.frame(method = AdWaldCI.df$method, x=AdWaldCI.df$x, LowerLimit = AdWaldCI.df$LATWx, UpperLimit = AdWaldCI.df$UATWx, LowerAbb = AdWaldCI.df$LABB, UpperAbb = AdWaldCI.df$UABB, ZWI = AdWaldCI.df$ZWI)

  Final.df= rbind(Generic.1,Generic.2,Generic.3,Generic.4,Generic.5, Generic.6)

  return(Final.df)
}
#####################################################################################
#' Plots the CI estimation of 6 adjusted methods (Wald, Wald-T, Likelihood, Score,
#' Logit-Wald, ArcSine)
#' @param x - Number of success
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @details  The plots of confidence intervals of using 6 adjusted methods
#' (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)
#' for \code{n} given \code{alp}, \code{x} and \code{h}
#' @family Adjusted methods of CI estimation  given x & n
#' @examples
#' x=5; n=5; alp=0.05; h=2
#' PlotciAAllx(x,n,alp,h)
#' @export
#10. Plot all methods
PlotciAAllx<-function(x,n,alp,h)
{
  if (missing(x)) stop("'x' is missing")
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (((class(x) != "integer") & (class(x) != "numeric")) || (x<0) || x>n || length(x)>1) stop("'x' has to be a positive integer between 0 and n")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h) >1|| h<0  || !(h%%1 ==0)) stop("'h' has to be an integer greater than or equal to 0")
  Abberation=ID=method=Value=LowerLimit=UpperLimit=LowerAbb=UpperAbb=ZWI=NULL

  ss1=ciAAllx(x,n,alp,h)
  id=1:nrow(ss1)
  ss= data.frame(ID=id,ss1)

  ll=subset(ss, LowerAbb=="YES")
  ul=subset(ss, UpperAbb=="YES")
  zl=subset(ss, ZWI=="YES")

  if (nrow(ll)>0) {
    ll=ll[,c(1,4)];
    ll$Abberation="Lower";
    colnames(ll)<-c("ID","Value","Abberation")}
  if (nrow(ul)>0){
    ul=ul[,c(1,5)]
    ul$Abberation="Upper"
    colnames(ul)<-c("ID","Value","Abberation")
  }
  if (nrow(zl)>0){
    zl=zl[,c(1,4)]
    zl$Abberation="ZWI"
    colnames(zl)<-c("ID","Value","Abberation")
  }
  ldf= rbind(ll,ul,zl)

  if(nrow(ldf)>0){
    oo= ggplot2::ggplot()+
      ggplot2::ggtitle("Confidence interval for adjusted methods") +
      ggplot2::labs(x = "Lower and Upper limits") +
      ggplot2::geom_errorbarh(data= ss,
                              ggplot2::aes(x = UpperLimit,y = ID,
                                           xmin = LowerLimit,
                                           xmax = UpperLimit,
                                           color= method),
                              size = 0.5)+
      ggplot2::geom_point(data=ldf,
                          ggplot2::aes(x=Value, y=ID,
                                       group = Abberation,shape=Abberation),
                          size = 4, fill = "red") +
      ggplot2::scale_fill_manual(values=c("blue", "cyan4", "red", "black", "orange","brown")) +
      ggplot2::scale_colour_manual(values=c("brown", "black", "blue", "cyan4", "red", "orange")) +
      ggplot2::scale_shape_manual(values=c(21,22,23))
  }
  else {
    oo=  ggplot2::ggplot()+
      ggplot2::ggtitle("Confidence interval for adjusted methods") +
      ggplot2::labs(x = "Lower and Upper limits") +
      ggplot2::geom_errorbarh(data= ss,
                              ggplot2::aes(x = UpperLimit,y = ID,
                                           xmin = LowerLimit,
                                           xmax = UpperLimit, color= method),
                              size = 0.5)
  }
  oo
}
#############################################
#' Plots the CI estimation of 6 adjusted methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine) grouped by x value
#' @param x - Number of success
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @details  The plots of confidence intervals of using 6 adjusted methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine) grouped by x for \code{n} given \code{alp}, \code{x} and \code{h}
#' @examples
#' x=5; n=5; alp=0.5;h=2
#' PlotciAAllxg(x,n,alp,h)
#' @export
#12.All methods plots with grouping
PlotciAAllxg<-function(x,n,alp,h)
{
  if (missing(x)) stop("'x' is missing")
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (((class(x) != "integer") & (class(x) != "numeric")) || (x<0) || x>n || length(x)>1) stop("'x' has to be a positive integer between 0 and n")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h) >1|| h<0  || !(h%%1 ==0)) stop("'h' has to be an integer greater than or equal to 0")
  Abberation=ID=method=Value=val1=val2=LowerLimit=UpperLimit=LowerAbb=UpperAbb=ZWI=NULL

  ss1=ciAAllx(x,n,alp,h)
  nss= ss1[order(ss1$x, (ss1$UpperLimit-ss1$LowerLimit)),]
  id=1:nrow(ss1)
  ss= data.frame(ID=id,nss)

  ll=subset(ss, LowerAbb=="YES")
  ul=subset(ss, UpperAbb=="YES")
  zl=subset(ss, ZWI=="YES")

  if (nrow(ll)>0) {
    ll=ll[,c(1,4)];
    ll$Abberation="Lower";
    colnames(ll)<-c("ID","Value","Abberation")}
  if (nrow(ul)>0){
    ul=ul[,c(1,5)]
    ul$Abberation="Upper"
    colnames(ul)<-c("ID","Value","Abberation")
  }
  if (nrow(zl)>0){
    zl=zl[,c(1,4)]
    zl$Abberation="ZWI"
    colnames(zl)<-c("ID","Value","Abberation")
  }
  ldf= rbind(ll,ul,zl)

  if((max(as.numeric(unique(ss$method)))-nrow(ss))==0){
    if(nrow(ldf)>0){
      oo= ggplot2::ggplot()+
        ggplot2::ggtitle("Confidence interval for adjusted methods sorted by x") +
        ggplot2::labs(x = "Lower and Upper limits") +
        ggplot2::geom_errorbarh(data= ss,
                                ggplot2::aes(x = UpperLimit,y = ID,
                                             xmin = LowerLimit,
                                             xmax = UpperLimit,
                                             color= method),
                                size = 0.5)+
        ggplot2::geom_point(data=ldf,
                            ggplot2::aes(x=Value, y=ID,
                                         group = Abberation,shape=Abberation),
                            size = 4, fill = "red") +
        ggplot2::scale_fill_manual(values=c("blue", "cyan4", "red", "black", "orange","brown")) +
        ggplot2::scale_colour_manual(values=c("brown", "black", "blue", "cyan4", "red", "orange")) +
        ggplot2::scale_shape_manual(values=c(21,22,23))
    }
    else {
      oo=  ggplot2::ggplot()+
        ggplot2::ggtitle("Confidence interval for adjusted methods sorted by x") +
        ggplot2::labs(x = "Lower and Upper limits") +
        ggplot2::geom_errorbarh(data= ss,
                                ggplot2::aes(x = UpperLimit,y = ID,
                                             xmin = LowerLimit,
                                             xmax = UpperLimit, color= method),
                                size = 0.5)
    }
    oo
  }
  else {

    ff= data.frame(val1=seq(0.5,max(ss$ID),by=(max(ss$ID)/(max(ss$x)+1))),val2=(0:max(ss$x)))

    if(nrow(ldf)>0){
      oo= ggplot2::ggplot()+
        ggplot2::ggtitle("Confidence interval for adjusted methods sorted by x") +
        ggplot2::labs(x = "Lower and Upper limits") +
        ggplot2::geom_errorbarh(data= ss,
                                ggplot2::aes(x = UpperLimit,y = ID,
                                             xmin = LowerLimit,
                                             xmax = UpperLimit,
                                             color= method),
                                size = 0.5)+
        ggplot2::geom_point(data=ldf,
                            ggplot2::aes(x=Value, y=ID,
                                         group = Abberation,shape=Abberation),
                            size = 4, fill = "red") +
        ggplot2::scale_fill_manual(values=c("blue", "cyan4", "red", "black", "orange","brown")) +
        ggplot2::scale_colour_manual(values=c("brown", "black", "blue", "cyan4", "red", "orange")) +
        ggplot2::scale_shape_manual(values=c(21,22,23))    +
        ggplot2::geom_hline(ggplot2::aes(yintercept=val1),data=ff) +
        ggplot2::geom_text(ggplot2::aes(0,val1,label = paste("x=", sep="", val2),hjust=1.1, vjust = -1), data=ff)
    }
    else {
      oo=  ggplot2::ggplot()+
        ggplot2::ggtitle("Confidence interval for adjusted methods given x & n sorted by x") +
        ggplot2::labs(x = "Lower and Upper limits") +
        ggplot2::geom_errorbarh(data= ss,
                                ggplot2::aes(x = UpperLimit,y = ID,
                                             xmin = LowerLimit,
                                             xmax = UpperLimit, color= method),
                                size = 0.5) +
        ggplot2::geom_hline(ggplot2::aes(yintercept=val1),data=ff) +
        ggplot2::geom_text(ggplot2::aes(0,val1,label = paste("x=", sep="", val2),hjust=1.1, vjust = -1), data=ff)
    }
    oo
  }
}
