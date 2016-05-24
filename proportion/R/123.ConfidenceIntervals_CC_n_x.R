#' Continuity corrected Wald method of CI estimation
#' @param x - Number of successes
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correction
#' @details  Wald-type interval (for all \eqn{x = 0, 1, 2 ..n}) using the test statistic
#' \eqn{(abs(phat-p)-c)/SE} where
#' \eqn{c > 0} is a constant for continuity correction
#' @return A dataframe with
#'  \item{x}{  Number of successes (positive samples)}
#'  \item{LCWx }{   CC-Wald Lower limit}
#'  \item{UCWx }{   CC-Wald Upper Limit}
#'  \item{LABB }{   CC-Wald Lower Abberation}
#'  \item{UABB }{   CC-Wald Upper Abberation}
#'  \item{ZWI }{   Zero Width Interval}
#' @family Continuity correction methods of CI estimation given x and n
#' @seealso \code{\link{prop.test} and \link{binom.test}} for equivalent base Stats R functionality,
#'    \code{\link[binom]{binom.confint}}  provides similar functionality for 11 methods,
#'    \code{\link[PropCIs]{wald2ci}} which provides multiple functions for CI calculation ,
#'    \code{\link[BlakerCI]{binom.blaker.limits}} which calculates Blaker CI which is not covered here and
#'    \code{\link[prevalence]{propCI}} which provides similar functionality.
#' @examples
#' x= 5; n=5; alp=0.05; c=1/(2*n)
#' ciCWDx(x,n,alp,c)
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
ciCWDx<-function(x,n,alp,c)
{
  if (missing(x)) stop("'x' is missing")
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if (((class(x) != "integer") & (class(x) != "numeric")) || (x<0) || x>n || length(x)>1) stop("'x' has to be a positive integer between 0 and n")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(c) != "integer") & (class(c) != "numeric") || length(c) >1 || c<0 ) stop("'c' has to be positive")

  ###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#WALD METHOD
pCWx=x/n
qCWx=1-pCWx
seCWx=sqrt(pCWx*qCWx/n)
LCWx=pCWx-((cv*seCWx)+c)
UCWx=pCWx+((cv*seCWx)+c)

if(LCWx<0) LABB="YES" else LABB="NO"
if(LCWx<0) LCWx=0

if(UCWx>1) UABB="YES" else UABB="NO"
if(UCWx>1) UCWx=1

if(UCWx-LCWx==0)ZWI="YES" else ZWI="NO"

return(data.frame(x,LCWx,UCWx,LABB,UABB,ZWI))
}
########################################################################################################
#' Continuity corrected Score method of CI estimation
#' @param x - Number of successes
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correction
#' @details  A score test approach using the
#' test statistic  \eqn{(abs(phat-p)-c)/SE}
#' where \eqn{c > 0} is a constant for continuity correction for
#' all \eqn{x = 0, 1, 2 ..n}
#' @return A dataframe with
#'  \item{x}{  Number of successes (positive samples)}
#'  \item{LCSx }{   Score Lower limit}
#'  \item{UCSx }{   Score Upper Limit}
#'  \item{LABB }{   Score Lower Abberation}
#'  \item{UABB }{   Score Upper Abberation}
#'  \item{ZWI }{   Zero Width Interval}
#' @family Continuity correction methods of CI estimation given x and n
#' @seealso \code{\link{prop.test} and \link{binom.test}} for equivalent base Stats R functionality,
#'    \code{\link[binom]{binom.confint}}  provides similar functionality for 11 methods,
#'    \code{\link[PropCIs]{wald2ci}} which provides multiple functions for CI calculation ,
#'    \code{\link[BlakerCI]{binom.blaker.limits}} which calculates Blaker CI which is not covered here and
#'    \code{\link[prevalence]{propCI}} which provides similar functionality.
#' @examples
#' x=5; n=5; alp=0.05; c=1/(2*n)
#' ciCSCx(x,n,alp,c)
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
ciCSCx<-function(x,n,alp,c)
{
  if (missing(x)) stop("'x' is missing")
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if (((class(x) != "integer") & (class(x) != "numeric")) || (x<0) || x>n || length(x)>1) stop("'x' has to be a positive integer between 0 and n")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>0.1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 0.1")
  if (c<=0 || c>(1/(2*n))) stop("'c' has to be positive and less than or equal to 1/(2*n)")


###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
cv1=(cv^2)/(2*n)
cv2= cv/(2*n)

#SCORE (WILSON) METHOD
pCSx=x/n
#qCSx=1-pCSx
seCS_Lx=sqrt((cv^2)-(4*n*(c+c^2))+(4*n*pCSx*(1-pCSx+(2*c))))	#Sq. root term of LL
seCS_Ux=sqrt((cv^2)+(4*n*(c-c^2))+(4*n*pCSx*(1-pCSx-(2*c))))	#Sq. root term of LL
LCSx=(n/(n+(cv)^2))*((pCSx-c+cv1)-(cv2*seCS_Lx))
UCSx=(n/(n+(cv)^2))*((pCSx+c+cv1)+(cv2*seCS_Ux))

if(LCSx<0) LABB="YES" else LABB="NO"
if(LCSx<0) LCSx=0

if(UCSx>1) UABB="YES" else UABB="NO"
if(UCSx>1) UCSx=1

if(UCSx-LCSx==0)ZWI="YES" else ZWI="NO"
return(data.frame(x,LCSx,UCSx,LABB,UABB,ZWI))
}
########################################################################################################
#' Continuity corrected ArcSine method of CI estimation
#' @param x - Number of successes
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correction
#' @details  Wald-type interval for the arcsine transformation using the test statistic
#' \eqn{(abs(sin^(-1)phat-sin^(-1)p)-c)/SE}
#'  where \eqn{c > 0} is a constant for continuity correction and for all \eqn{x = 0, 1, 2 ..n}
#' @return A dataframe with
#'  \item{x}{  Number of successes (positive samples)}
#'  \item{LCAx }{   ArcSine Lower limit}
#'  \item{UCAx }{   ArcSine Upper Limit}
#'  \item{LABB }{   ArcSine Lower Abberation}
#'  \item{UABB }{   ArcSine Upper Abberation}
#'  \item{ZWI }{   Zero Width Interval}
#' @family Continuity correction methods given x and n
#' @seealso \code{\link{prop.test} and \link{binom.test}} for equivalent base Stats R functionality,
#'    \code{\link[binom]{binom.confint}}  provides similar functionality for 11 methods,
#'    \code{\link[PropCIs]{wald2ci}} which provides multiple functions for CI calculation ,
#'    \code{\link[BlakerCI]{binom.blaker.limits}} which calculates Blaker CI which is not covered here and
#'    \code{\link[prevalence]{propCI}} which provides similar functionality.
#' @examples
#' x=5; n=5; alp=0.05;c=1/2*n
#' ciCASx(x,n,alp,c)
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
#3.ARC-SINE
ciCASx<-function(x,n,alp,c)
{
  if (missing(x)) stop("'x' is missing")
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if (((class(x) != "integer") & (class(x) != "numeric")) || (x<0) || x>n || length(x)>1) stop("'x' has to be a positive integer between 0 and n")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(c) != "integer") & (class(c) != "numeric") || length(c) >1 || c<0 ) stop("'c' has to be positive")

###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#ARC-SINE METHOD
pCAx=x/n
#qCAx=1-pCAx
seCAx=cv/sqrt(4*n)
LCAx=(sin(asin(sqrt(pCAx))-seCAx-c))^2
UCAx=(sin(asin(sqrt(pCAx))+seCAx+c))^2

if(LCAx<0) LABB ="YES" else LABB ="NO"
if(LCAx<0) LCAx=0

if(UCAx>1) UABB ="YES" else UABB="NO"
if(UCAx>1) UCAx=1

if(UCAx-LCAx==0)ZWI="YES" else ZWI="NO"

return(data.frame(x,LCAx,UCAx,LABB,UABB,ZWI))
}
########################################################################################################
#' Continuity corrected Logit-Wald method of CI estimation
#' @param x - Number of successes
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correction
#' @details  Wald-type interval for the logit transformation of the parameter \code{p}
#' using the test statistic
#' \eqn{(abs(L(phat)-L(p))-c)/SE}
#' where \eqn{c > 0} is a constant for continuity correction and \eqn{L(y) = log(y/1-y)}
#' for all \eqn{x = 0, 1, 2 ..n}. Boundary modifications when \eqn{x = 0} or \eqn{x = n}
#' using Exact method values.
#' @return A dataframe with
#'  \item{x}{  Number of successes (positive samples)}
#'  \item{LCLTx }{   Logit Wald Lower limit}
#'  \item{UCLTx }{   Logit Wald Upper Limit}
#'  \item{LABB }{   Logit Wald Lower Abberation}
#'  \item{UABB }{   Logit Wald Upper Abberation}
#'  \item{ZWI }{   Zero Width Interval}
#' @family Continuity correction methods of CI estimation given x and n
#' @seealso \code{\link{prop.test} and \link{binom.test}} for equivalent base Stats R functionality,
#'    \code{\link[binom]{binom.confint}}  provides similar functionality for 11 methods,
#'    \code{\link[PropCIs]{wald2ci}} which provides multiple functions for CI calculation ,
#'    \code{\link[BlakerCI]{binom.blaker.limits}} which calculates Blaker CI which is not covered here and
#'    \code{\link[prevalence]{propCI}} which provides similar functionality.
#' @examples
#' x=5; n=5; alp=0.05;c=1/2*n
#' ciCLTx(x,n,alp,c)
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
#4.LOGIT-WALD
ciCLTx<-function(x,n,alp,c)
{
  if (missing(x)) stop("'x' is missing")
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if (((class(x) != "integer") & (class(x) != "numeric")) || (x<0) || x>n || length(x)>1) stop("'x' has to be a positive integer between 0 and n")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(c) != "integer") & (class(c) != "numeric") || length(c) >1 || c<0 ) stop("'c' has to be positive")

###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#LOGIT-WALD METHOD
if(x==0)
{
pCLTx=0
qCLTx=1
LCLTx= 0
UCLTx = 1-((alp/2)^(1/n))
}
if(x==n)
{
pCLTx=1
qCLTx=0
LCLTx= (alp/2)^(1/n)
UCLTx=1
}
if(0 < x && x < n)
{
lgiti=function(t) exp(t)/(1+exp(t))	#LOGIT INVERSE
pCLTx=x/n
qCLTx=1-pCLTx
lgitx=log(pCLTx/qCLTx)
seCLTx=sqrt(pCLTx*qCLTx*n)
LCLTx=lgiti(lgitx-(cv/seCLTx)-c)
UCLTx=lgiti(lgitx+(cv/seCLTx)+c)
}

if(LCLTx<0) LABB="YES" else LABB="NO"
if(LCLTx<0) LCLTx=0

if(UCLTx>1) UABB="YES" else UABB="NO"
if(UCLTx>1) UCLTx=1

if(UCLTx-LCLTx==0)ZWI="YES" else ZWI="NO"

return(data.frame(x,LCLTx,UCLTx,LABB,UABB,ZWI))
}
############################################################################
#' Continuity corrected Wald-T method of CI estimation
#' @param x - Number of successes
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correction
#' @details  Approximate method based on a t_approximation of the standardized point estimator
#' using the test statistic
#' \eqn{(abs(phat-p)-c)/SE}
#' where \eqn{c > 0} is a constant for continuity correction for all \eqn{x = 0, 1, 2 ..n}.
#' Boundary modifications when \eqn{x = 0} or \eqn{x = n} using Wald adjustment method with
#' \eqn{h = 2}.
#' @return A dataframe with
#'  \item{x}{  Number of successes (positive samples)}
#'  \item{LCTWx }{   T-Wald Lower limit}
#'  \item{UCTWx }{   T-Wald Upper Limit}
#'  \item{LABB }{   T-Wald Lower Abberation}
#'  \item{UABB }{   T-Wald Upper Abberation}
#'  \item{ZWI }{   Zero Width Interval}
#' @family Continuity correction methods of CI estimation given x and n
#' @seealso \code{\link{prop.test} and \link{binom.test}} for equivalent base Stats R functionality,
#'    \code{\link[binom]{binom.confint}}  provides similar functionality for 11 methods,
#'    \code{\link[PropCIs]{wald2ci}} which provides multiple functions for CI calculation ,
#'    \code{\link[BlakerCI]{binom.blaker.limits}} which calculates Blaker CI which is not covered here and
#'    \code{\link[prevalence]{propCI}} which provides similar functionality.
#' @examples
#' x=5; n=5; alp=0.05;c=1/2*n
#' ciCTWx(x,n,alp,c)
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
#5.T-WALD
ciCTWx<-function(x,n,alp,c)
{
  if (missing(x)) stop("'x' is missing")
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if (((class(x) != "integer") & (class(x) != "numeric")) || (x<0) || x>n || length(x)>1) stop("'x' has to be a positive integer between 0 and n")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(c) != "integer") & (class(c) != "numeric") || length(c) >1 || c<0 ) stop("'c' has to be positive")

#MODIFIED_t-WALD METHOD
if(x==0||x==n)
{
pCTWx=(x+2)/(n+4)
#qCTWx=1-pCTWx
}else
{
pCTWx=x/n
#qCTWx=1-pCTWx
}
f1=function(p,n) p*(1-p)/n
f2=function(p,n) (p*(1-p)/(n^3))+(p+((6*n)-7)*(p^2)+(4*(n-1)*(n-3)*(p^3))-(2*(n-1)*((2*n)-3)*(p^4)))/(n^5)-(2*(p+((2*n)-3)*(p^2)-2*(n-1)*(p^3)))/(n^4)
DOFx=2*((f1(pCTWx,n))^2)/f2(pCTWx,n)
cvx=qt(1-(alp/2), df=DOFx)
seCTWx=cvx*sqrt(f1(pCTWx,n))
LCTWx=pCTWx-(seCTWx+c)
UCTWx=pCTWx+(seCTWx+c)

if(LCTWx<0) LABB="YES" else LABB="NO"
if(LCTWx<0) LCTWx=0

if(UCTWx>1) UABB="YES" else  UABB="NO"
if(UCTWx>1) UCTWx=1

if(UCTWx-LCTWx==0)ZWI ="YES" else ZWI="NO"

return(data.frame(x,LCTWx,UCTWx,LABB,UABB,ZWI))
}
####################################################################################
#' CI estimation of 5 continuity corrected methods (Wald, Wald-T, Score, Logit-Wald, ArcSine)
#' @param x - Number of sucess
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correction
#' @details  The Confidence Interval of 5 continuity corrected methods (Wald, Wald-T, Score, Logit-Wald, ArcSine) for \code{n} given \code{alp} and \code{x}
#' @return A dataframe with
#'  \item{method }{- Name of the method}
#'  \item{x }{- Number of successes (positive samples)}
#'  \item{LLT }{ - Lower limit}
#'  \item{ULT }{ - Upper Limit}
#'  \item{LABB }{ - Lower Abberation}
#'  \item{UABB }{ - Upper Abberation}
#'  \item{ZWI }{ - Zero Width Interval}
#' @family Continuity correction methods of CI estimation given x and n
#' @seealso \code{\link{prop.test} and \link{binom.test}} for equivalent base Stats R functionality,
#'    \code{\link[binom]{binom.confint}}  provides similar functionality for 11 methods,
#'    \code{\link[PropCIs]{wald2ci}} which provides multiple functions for CI calculation ,
#'    \code{\link[BlakerCI]{binom.blaker.limits}} which calculates Blaker CI which is not covered here and
#'    \code{\link[prevalence]{propCI}} which provides similar functionality.
#' @examples
#' x=5; n=5; alp=0.05;c=1/(2*n)
#' ciCAllx(x,n,alp,c)
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
#6.All methods
ciCAllx<-function(x,n,alp,c)
{
  if (missing(x)) stop("'x' is missing")
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if (((class(x) != "integer") & (class(x) != "numeric")) || (x<0) || x>n || length(x)>1) stop("'x' has to be a positive integer between 0 and n")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>.1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and .1")
  if ((class(c) != "integer") & (class(c) != "numeric") || length(c) >1 || c<0 ) stop("'c' has to be positive")

  #### Calling functions and creating df
  WaldCI.df    = ciCWDx(x,n,alp,c)
  ArcSineCI.df = ciCASx(x,n,alp,c)
  ScoreCI.df   = ciCSCx(x,n,alp,c)
  WaldLCI.df   = ciCLTx(x,n,alp,c)
  WaldTCI.df   = ciCTWx(x,n,alp,c)

  WaldCI.df$method    = as.factor("Wald")
  ArcSineCI.df$method = as.factor("ArcSine")
  WaldLCI.df$method   = as.factor("Logit Wald")
  ScoreCI.df$method   = as.factor("Score")
  WaldTCI.df$method   = as.factor("Wald-T")

  Generic.1 = data.frame(method = WaldCI.df$method, x=WaldCI.df$x, LowerLimit = WaldCI.df$LCWx, UpperLimit = WaldCI.df$UCWx, LowerAbb = WaldCI.df$LABB, UpperAbb = WaldCI.df$UABB, ZWI = WaldCI.df$ZWI)
  Generic.2 = data.frame(method = ArcSineCI.df$method, x=ArcSineCI.df$x, LowerLimit = ArcSineCI.df$LCAx, UpperLimit = ArcSineCI.df$UCAx, LowerAbb = ArcSineCI.df$LABB, UpperAbb = ArcSineCI.df$UABB, ZWI = ArcSineCI.df$ZWI)
  Generic.4 = data.frame(method = ScoreCI.df$method, x=ScoreCI.df$x, LowerLimit = ScoreCI.df$LCSx, UpperLimit = ScoreCI.df$UCSx, LowerAbb = ScoreCI.df$LABB, UpperAbb = ScoreCI.df$UABB, ZWI = ScoreCI.df$ZWI)
  Generic.5 = data.frame(method = WaldLCI.df$method, x=WaldLCI.df$x, LowerLimit = WaldLCI.df$LCLTx, UpperLimit = WaldLCI.df$UCLTx, LowerAbb = WaldLCI.df$LABB, UpperAbb = WaldLCI.df$UABB, ZWI = WaldLCI.df$ZWI)
  Generic.6 = data.frame(method = WaldTCI.df$method, x=WaldTCI.df$x, LowerLimit = WaldTCI.df$LCTWx, UpperLimit = WaldTCI.df$UCTWx, LowerAbb = WaldTCI.df$LABB, UpperAbb = WaldTCI.df$UABB, ZWI = WaldTCI.df$ZWI)

  Final.df= rbind(Generic.1,Generic.2,Generic.4,Generic.5, Generic.6)

  return(Final.df)
}
#####################################################################################
#' Plots the CI estimation of 5 continuity corrected methods (Wald, Wald-T, Score, Logit-Wald, ArcSine) given n, alp and x with continuity correction c
#' @param x - Number of sucess
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correction
#' @details  Plots the Confidence Interval of 5 continuity corrected methods (Wald, Wald-T, Score, Logit-Wald, ArcSine) for \code{n} given \code{alp}, \code{x} and and \code{c}
#' @family Continuity correction methods of CI estimation given x and n
#' @examples
#' x=5; n=5; alp=0.05;c=1/(2*n)
#' PlotciCAllx(x,n,alp,c)
#' @export
#7. Plot all methods
PlotciCAllx<-function(x,n,alp,c)
{
  if (missing(x)) stop("'x' is missing")
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if (((class(x) != "integer") & (class(x) != "numeric")) || (x<0) || x>n || length(x)>1) stop("'x' has to be a positive integer between 0 and n")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(c) != "integer") & (class(c) != "numeric") || length(c) >1 || c<0 ) stop("'c' has to be positive")
  Abberation=ID=method=Value=LowerLimit=UpperLimit=LowerAbb=UpperAbb=ZWI=NULL

  ss1=ciCAllx(x,n,alp,c)
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
      ggplot2::ggtitle("Confidence interval for continuity correction methods") +
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
      ggplot2::ggtitle("Confidence interval for continuity correction methods") +
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
#' Plots the CI estimation of 5 continuity corrected methods (Wald, Wald-T, Score, Logit-Wald, ArcSine) given x & n grouped by x value
#' @param x - Number of sucess
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correction
#' @details  Plots he Confidence Interval of 5 continuity corrected methods (Wald, Wald-T, Score, Logit-Wald, ArcSine) for \code{n} given \code{alp}, \code{c} and \code{x} grouped by x
#' @family Continuity correction methods of CI estimation given x and n
#' @examples
#' x=5; n=5; alp=0.05;c=1/(2*n)
#' PlotciCAllxg(x,n,alp,c)
#' @export
#8.All methods plots with grouping for 5 continuity corrected methods (Wald, Wald-T, Score, Logit-Wald, ArcSine)
PlotciCAllxg<-function(x,n,alp,c)
{
  if (missing(x)) stop("'x' is missing")
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if (((class(x) != "integer") & (class(x) != "numeric")) || (x<0) || x>n || length(x)>1) stop("'x' has to be a positive integer between 0 and n")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(c) != "integer") & (class(c) != "numeric") || length(c) >1 || c<0 ) stop("'c' has to be positive")
  Abberation=ID=method=Value=val1=val2=LowerLimit=UpperLimit=LowerAbb=UpperAbb=ZWI=NULL

  ss1=ciCAllx(x,n,alp,c)
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
        ggplot2::ggtitle("Confidence interval for continuity correction methods given x & n sorted by x") +
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
        ggplot2::ggtitle("Confidence interval for continuity correction methods given x & n sorted by x") +
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
        ggplot2::ggtitle("Confidence interval for continuity correction methods given x & n sorted by x") +
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
        ggplot2::ggtitle("Confidence interval for continuity correction methods given x & n sorted by x") +
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
