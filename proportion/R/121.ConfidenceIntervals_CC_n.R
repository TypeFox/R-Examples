#' Continuity corrected Wald method of CI estimation
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correction
#' @details  Wald-type interval (for all \eqn{x = 0, 1, 2 ..n}) using the test statistic \eqn{(abs(phat-p)-c)/SE} where
#' \eqn{c > 0} is a constant for continuity correction
#' @return A dataframe with
#'  \item{x}{  Number of successes (positive samples)}
#'  \item{LCW }{   Wald Lower limit}
#'  \item{UCW }{   Wald Upper Limit}
#'  \item{LABB }{   Wald Lower Abberation}
#'  \item{UABB }{   Wald Upper Abberation}
#'  \item{ZWI }{   Zero Width Interval}
#' @family Continuity correction methods of CI estimation
#' @seealso \code{\link{prop.test} and \link{binom.test}} for equivalent base Stats R functionality,
#'    \code{\link[binom]{binom.confint}}  provides similar functionality for 11 methods,
#'    \code{\link[PropCIs]{wald2ci}} which provides multiple functions for CI calculation ,
#'    \code{\link[BlakerCI]{binom.blaker.limits}} which calculates Blaker CI which is not covered here and
#'    \code{\link[prevalence]{propCI}} which provides similar functionality.
#' @examples
#' n=5; alp=0.05;c=1/(2*n)
#' ciCWD(n,alp,c)
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
ciCWD<-function(n,alp,c)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if ((class(c) != "integer") & (class(c) != "numeric") || length(c) >1 || c<0 ) stop("'c' has to be positive")

####INPUT n
x=0:n
k=n+1
####INITIALIZATIONS
pCW=0
qCW=0
seCW=0
LCW=0
UCW=0
LABB=0
UABB=0
ZWI=0
###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#WALD METHOD
for(i in 1:k)
{
pCW[i]=x[i]/n
qCW[i]=1-pCW[i]
seCW[i]=sqrt(pCW[i]*qCW[i]/n)
LCW[i]=pCW[i]-((cv*seCW[i])+c)
UCW[i]=pCW[i]+((cv*seCW[i])+c)

if(LCW[i]<0) LABB[i]="YES" else LABB[i]="NO"
if(LCW[i]<0) LCW[i]=0

if(UCW[i]>1) UABB[i]="YES" else UABB[i]="NO"
if(UCW[i]>1) UCW[i]=1

if(UCW[i]-LCW[i]==0)ZWI[i]="YES" else ZWI[i]="NO"
}
return(data.frame(x,LCW,UCW,LABB,UABB,ZWI))
}
########################################################################################################
##############################################################################
#' Continuity corrected Score method of CI estimation
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correction
#' @details  A score test approach using the
#' test statistic  \eqn{(abs(phat-p)-c)/SE}
#' where \eqn{c > 0} is a constant for continuity correction for
#' all \eqn{x = 0, 1, 2 ..n}
#' @return A dataframe with
#'  \item{x}{  Number of successes (positive samples)}
#'  \item{LCS }{   Score Lower limit}
#'  \item{UCS }{   Score Upper Limit}
#'  \item{LABB }{   Score Lower Abberation}
#'  \item{UABB }{   Score Upper Abberation}
#'  \item{ZWI }{   Zero Width Interval}
#' @family Continuity correction methods of CI estimation
#' @seealso \code{\link{prop.test} and \link{binom.test}} for equivalent base Stats R functionality,
#'    \code{\link[binom]{binom.confint}}  provides similar functionality for 11 methods,
#'    \code{\link[PropCIs]{wald2ci}} which provides multiple functions for CI calculation ,
#'    \code{\link[BlakerCI]{binom.blaker.limits}} which calculates Blaker CI which is not covered here and
#'    \code{\link[prevalence]{propCI}} which provides similar functionality.
#' @examples
#' n=5; alp=0.05;c=1/(2*n)
#' ciCSC(n,alp,c)
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
ciCSC<-function(n,alp,c)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (c<=0 || c>(1/(2*n))) stop("'c' has to be positive and less than or equal to 1/(2*n)")

  ####INPUT n
  x=0:n
  k=n+1
  ####INITIALIZATIONS
  pCS=0
  qCS=0
  seCS_L=0
  seCS_U=0
  LCS=0
  UCS=0
  LABB=0
  UABB=0
  ZWI=0

  ###CRITICAL VALUES
  cv=qnorm(1-(alp/2), mean = 0, sd = 1)
  cv1=(cv^2)/(2*n)
  cv2= cv/(2*n)

  #SCORE (WILSON) METHOD
  for(i in 1:k)
  {
    pCS[i]=x[i]/n
    qCS[i]=1-pCS[i]
    seCS_L[i]=sqrt((cv^2)-(4*n*(c+c^2))+(4*n*pCS[i]*(1-pCS[i]+(2*c))))	#Sq. root term of LL
    seCS_U[i]=sqrt((cv^2)+(4*n*(c-c^2))+(4*n*pCS[i]*(1-pCS[i]-(2*c))))	#Sq. root term of LL
    LCS[i]=(n/(n+(cv)^2))*((pCS[i]-c+cv1)-(cv2*seCS_L[i]))
    UCS[i]=(n/(n+(cv)^2))*((pCS[i]+c+cv1)+(cv2*seCS_U[i]))

    if(LCS[i]<0) LABB[i]="YES" else LABB[i]="NO"
    if(LCS[i]<0) LCS[i]=0

    if(UCS[i]>1) UABB[i]="YES" else UABB[i]="NO"
    if(UCS[i]>1) UCS[i]=1

    if(UCS[i]-LCS[i]==0)ZWI[i]="YES" else ZWI[i]="NO"
  }
  return(data.frame(x,LCS,UCS,LABB,UABB,ZWI))
}
########################################################################################################
#' Continuity corrected ArcSine method of CI estimation
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correction
#' @details  Wald-type interval for the arcsine transformation using the test statistic
#' \eqn{(abs(sin^(-1)phat-sin^(-1)p)-c)/SE}
#'  where \eqn{c > 0} is a constant for continuity correction and for all \eqn{x = 0, 1, 2 ..n}
#' @return A dataframe with
#'  \item{x}{  Number of successes (positive samples)}
#'  \item{LCA }{   ArcSine Lower limit}
#'  \item{UCA }{   ArcSine Upper Limit}
#'  \item{LABB }{   ArcSine Lower Abberation}
#'  \item{UABB }{   ArcSine Upper Abberation}
#'  \item{ZWI }{   Zero Width Interval}
#' @family Continuity correction methods of CI estimation
#' @seealso \code{\link{prop.test} and \link{binom.test}} for equivalent base Stats R functionality,
#'    \code{\link[binom]{binom.confint}}  provides similar functionality for 11 methods,
#'    \code{\link[PropCIs]{wald2ci}} which provides multiple functions for CI calculation ,
#'    \code{\link[BlakerCI]{binom.blaker.limits}} which calculates Blaker CI which is not covered here and
#'    \code{\link[prevalence]{propCI}} which provides similar functionality.
#' @examples
#' n=5; alp=0.05;c=1/(2*n)
#' ciCAS(n,alp,c)
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
ciCAS<-function(n,alp,c)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if ((class(c) != "integer") & (class(c) != "numeric") || length(c) >1 || c<0 ) stop("'c' has to be positive")

####INPUT n
x=0:n
k=n+1
####INITIALIZATIONS
pCA=0
qCA=0
seCA=0
LCA=0
UCA=0
LABB=0
UABB=0
ZWI=0

###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#ARC-SINE METHOD
for(i in 1:k)
{
pCA[i]=x[i]/n
qCA[i]=1-pCA[i]
seCA[i]=cv/sqrt(4*n)
LCA[i]=(sin(asin(sqrt(pCA[i]))-seCA[i]-c))^2
UCA[i]=(sin(asin(sqrt(pCA[i]))+seCA[i]+c))^2

if(LCA[i]<0) LABB[i]="YES" else LABB[i]="NO"
if(LCA[i]<0) LCA[i]=0

if(UCA[i]>1) UABB[i]="YES" else UABB[i]="NO"
if(UCA[i]>1) UCA[i]=1

if(UCA[i]-LCA[i]==0)ZWI[i]="YES" else ZWI[i]="NO"
}
return(data.frame(x,LCA,UCA,LABB,UABB,ZWI))
}

########################################################################################################
#' Continuity corrected Logit Wald method of CI estimation
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
#'  \item{LCLT }{   Logit Wald Lower limit}
#'  \item{UCLT }{   Logit Wald Upper Limit}
#'  \item{LABB }{   Logit Wald Lower Abberation}
#'  \item{UABB }{   Logit Wald Upper Abberation}
#'  \item{ZWI }{   Zero Width Interval}
#' @family Continuity correction methods of CI estimation
#' @seealso \code{\link{prop.test} and \link{binom.test}} for equivalent base Stats R functionality,
#'    \code{\link[binom]{binom.confint}}  provides similar functionality for 11 methods,
#'    \code{\link[PropCIs]{wald2ci}} which provides multiple functions for CI calculation ,
#'    \code{\link[BlakerCI]{binom.blaker.limits}} which calculates Blaker CI which is not covered here and
#'    \code{\link[prevalence]{propCI}} which provides similar functionality.
#' @examples
#' n=5; alp=0.05;c=1/(2*n)
#' ciCLT(n,alp,c)
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
ciCLT<-function(n,alp,c)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if ((class(c) != "integer") & (class(c) != "numeric") || length(c) >1 || c<0 ) stop("'c' has to be positive")

####INPUT n
x=0:n
k=n+1
####INITIALIZATIONS
pCLT=0
qCLT=0
seCLT=0
lgit=0
LCLT=0
UCLT=0
LABB=0
UABB=0
ZWI=0
###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#LOGIT-WALD METHOD
pCLT[1]=0
qCLT[1]=1
LCLT[1] = 0
UCLT[1] = 1-((alp/2)^(1/n))

pCLT[k]=1
qCLT[k]=0
LCLT[k]= (alp/2)^(1/n)
UCLT[k]=1

lgiti=function(t) exp(t)/(1+exp(t))	#LOGIT INVERSE
for(j in 1:(k-2))
{
pCLT[j+1]=x[j+1]/n
qCLT[j+1]=1-pCLT[j+1]
lgit[j+1]=log(pCLT[j+1]/qCLT[j+1])
seCLT[j+1]=sqrt(pCLT[j+1]*qCLT[j+1]*n)
LCLT[j+1]=lgiti(lgit[j+1]-(cv/seCLT[j+1])-c)
UCLT[j+1]=lgiti(lgit[j+1]+(cv/seCLT[j+1])+c)
}
for(i in 1:k)
{
if(LCLT[i]<0) LABB[i]="YES" else LABB[i]="NO"
if(LCLT[i]<0) LCLT[i]=0

if(UCLT[i]>1) UABB[i]="YES" else UABB[i]="NO"
if(UCLT[i]>1) UCLT[i]=1

if(UCLT[i]-LCLT[i]==0)ZWI[i]="YES" else ZWI[i]="NO"
}
return(data.frame(x,LCLT,UCLT,LABB,UABB,ZWI))
}
############################################################################
#' Continuity corrected  Wald-T method of CI estimation
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correction
#' @details  Approximate method based on a t_approximation of the standardized
#' point estimator using the test statistic
#' \eqn{(abs(phat-p)-c)/SE}
#' where \eqn{c > 0} is a constant for continuity correction for all \eqn{x = 0, 1, 2 ..n}.
#' Boundary modifications when \eqn{x = 0} or \eqn{x = n} using Wald adjustment method with
#' \eqn{h = 2}.
#' @return A dataframe with
#'  \item{x}{  Number of successes (positive samples)}
#'  \item{LCTW }{   T-Wald Lower limit}
#'  \item{UCTW }{   T-Wald Upper Limit}
#'  \item{LABB }{   T-Wald Lower Abberation}
#'  \item{UABB }{   T-Wald Upper Abberation}
#'  \item{ZWI }{   Zero Width Interval}
#' @family Continuity correction methods of CI estimation
#' @seealso \code{\link{prop.test} and \link{binom.test}} for equivalent base Stats R functionality,
#'    \code{\link[binom]{binom.confint}}  provides similar functionality for 11 methods,
#'    \code{\link[PropCIs]{wald2ci}} which provides multiple functions for CI calculation ,
#'    \code{\link[BlakerCI]{binom.blaker.limits}} which calculates Blaker CI which is not covered here and
#'    \code{\link[prevalence]{propCI}} which provides similar functionality.
#' @examples
#' n=5; alp=0.05;c=1/(2*n)
#' ciCTW(n,alp,c)
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
ciCTW<-function(n,alp,c)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if ((class(c) != "integer") & (class(c) != "numeric") || length(c) >1 || c<0 ) stop("'c' has to be positive")


####INPUT n
x=0:n
k=n+1
####INITIALIZATIONS
pCTW=0
qCTW=0
seCTW=0
LCTW=0
UCTW=0
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
pCTW[i]=(x[i]+2)/(n+4)
qCTW[i]=1-pCTW[i]
}else
{
pCTW[i]=x[i]/n
qCTW[i]=1-pCTW[i]
}
f1=function(p,n) p*(1-p)/n
f2=function(p,n) (p*(1-p)/(n^3))+(p+((6*n)-7)*(p^2)+(4*(n-1)*(n-3)*(p^3))-(2*(n-1)*((2*n)-3)*(p^4)))/(n^5)-(2*(p+((2*n)-3)*(p^2)-2*(n-1)*(p^3)))/(n^4)
DOF[i]=2*((f1(pCTW[i],n))^2)/f2(pCTW[i],n)
cv[i]=qt(1-(alp/2), df=DOF[i])
seCTW[i]=cv[i]*sqrt(f1(pCTW[i],n))
LCTW[i]=pCTW[i]-(seCTW[i]+c)
UCTW[i]=pCTW[i]+(seCTW[i]+c)

if(LCTW[i]<0) LABB[i]="YES" else LABB[i]="NO"
if(LCTW[i]<0) LCTW[i]=0

if(UCTW[i]>1) UABB[i]="YES" else  UABB[i]="NO"
if(UCTW[i]>1) UCTW[i]=1

if(UCTW[i]-LCTW[i]==0)ZWI[i]="YES" else ZWI[i]="NO"
}
return(data.frame(x,LCTW,UCTW,LABB,UABB,ZWI))
}
#####################################################################
#' CI estimation of 5 continuity corrected methods (Wald, Wald-T, Score, Logit-Wald, ArcSine)
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correction
#' @details  The Confidence Interval on 5 continuity corrected methods (Wald, Wald-T, Score, Logit-Wald, ArcSine) for \code{n} given \code{alp} along with Continuity correction \code{c}
#' @return A dataframe with
##'  \item{method }{- Name of the method}
##'  \item{x }{- Number of successes (positive samples)}
##'  \item{LLT }{ - Lower limit}
##'  \item{ULT }{ - Upper Limit}
##'  \item{LABB }{ - Lower Abberation}
##'  \item{UABB }{ - Upper Abberation}
##'  \item{ZWI }{ - Zero Width Interval}
#' @family Continuity correction methods of CI estimation
#' @seealso \code{\link{prop.test} and \link{binom.test}} for equivalent base Stats R functionality,
#'    \code{\link[binom]{binom.confint}}  provides similar functionality for 11 methods,
#'    \code{\link[PropCIs]{wald2ci}} which provides multiple functions for CI calculation ,
#'    \code{\link[BlakerCI]{binom.blaker.limits}} which calculates Blaker CI which is not covered here and
#'    \code{\link[prevalence]{propCI}} which provides similar functionality.
#' @examples
#' n=5; alp=0.05;c=1/(2*n)
#' ciCAll(n,alp,c)
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
ciCAll<-function(n,alp,c)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (c<=0 || c>(1/(2*n))) stop("'c' has to be positive and less than or equal to 1/(2*n)")

#### Calling functions and creating df
  WaldCI.df    = ciCWD(n,alp,c)
  ArcSineCI.df = ciCAS(n,alp,c)
  ScoreCI.df   = ciCSC(n,alp,c)
  WaldLCI.df   = ciCLT(n,alp,c)
  WaldTCI.df   = ciCTW(n,alp,c)

  WaldCI.df$method    = as.factor("CC-Wald")
  ArcSineCI.df$method = as.factor("CC-ArcSine")
  WaldLCI.df$method   = as.factor("CC-Logit-Wald")
  ScoreCI.df$method   = as.factor("CC-Score")
  WaldTCI.df$method   = as.factor("CC-Wald-T")

  Generic.1 = data.frame(method = WaldCI.df$method, x=WaldCI.df$x, LowerLimit = WaldCI.df$LCW, UpperLimit = WaldCI.df$UCW, LowerAbb = WaldCI.df$LABB, UpperAbb = WaldCI.df$UABB, ZWI = WaldCI.df$ZWI)
  Generic.2 = data.frame(method = ArcSineCI.df$method, x=ArcSineCI.df$x, LowerLimit = ArcSineCI.df$LCA, UpperLimit = ArcSineCI.df$UCA, LowerAbb = ArcSineCI.df$LABB, UpperAbb = ArcSineCI.df$UABB, ZWI = ArcSineCI.df$ZWI)
  Generic.4 = data.frame(method = ScoreCI.df$method, x=ScoreCI.df$x, LowerLimit = ScoreCI.df$LCS, UpperLimit = ScoreCI.df$UCS, LowerAbb = ScoreCI.df$LABB, UpperAbb = ScoreCI.df$UABB, ZWI = ScoreCI.df$ZWI)
  Generic.5 = data.frame(method = WaldLCI.df$method, x=WaldLCI.df$x, LowerLimit = WaldLCI.df$LCLT, UpperLimit = WaldLCI.df$UCLT, LowerAbb = WaldLCI.df$LABB, UpperAbb = WaldLCI.df$UABB, ZWI = WaldLCI.df$ZWI)
  Generic.6 = data.frame(method = WaldTCI.df$method, x=WaldTCI.df$x, LowerLimit = WaldTCI.df$LCTW, UpperLimit = WaldTCI.df$UCTW, LowerAbb = WaldTCI.df$LABB, UpperAbb = WaldTCI.df$UABB, ZWI = WaldTCI.df$ZWI)

  Final.df= rbind(Generic.1,Generic.2,Generic.4,Generic.5, Generic.6)

  return(Final.df)
}
