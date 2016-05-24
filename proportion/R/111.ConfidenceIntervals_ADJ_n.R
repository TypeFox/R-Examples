#' Adjusted Wald method of CI estimation
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - adding factor
#' @details  Given data \code{x} and \code{n} are modified as \eqn{x + h} and \eqn{n + (2*h)}
#' respectively, where \eqn{h > 0} then Wald-type interval is applied for all \eqn{x = 0, 1, 2 ..n.}
#' @return A dataframe with
#'  \item{x}{  Number of successes (positive samples)}
#'  \item{LAWD }{   Wald Lower limit}
#'  \item{UAWD }{   Wald Upper Limit}
#'  \item{LABB }{   Wald Lower Abberation}
#'  \item{UABB }{   Wald Upper Abberation}
#'  \item{ZWI }{   Zero Width Interval}
#' @family Adjusted methods of CI estimation
#' @seealso \code{\link{prop.test} and \link{binom.test}} for equivalent base Stats R functionality,
#'    \code{\link[binom]{binom.confint}}  provides similar functionality for 11 methods,
#'    \code{\link[PropCIs]{wald2ci}} which provides multiple functions for CI calculation ,
#'    \code{\link[BlakerCI]{binom.blaker.limits}} which calculates Blaker CI which is not covered here and
#'    \code{\link[prevalence]{propCI}} which provides similar functionality.
#' @examples
#' n=5; alp=0.05;h=2
#' ciAWD(n,alp,h)
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
ciAWD<-function(n,alp,h)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h) >1|| h<0  ) stop("'h' has to be greater than or equal to 0")

####INPUT n
 x=0:n
 k=n+1
 y=x+h
 n1=n+(2*h)
####INITIALIZATIONS
pAW=0
qAW=0
seAW=0
LAWD=0
UAWD=0
LABB=0
UABB=0
ZWI=0

###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#WALD METHOD
for(i in 1:k)
{
pAW[i]=y[i]/n1
qAW[i]=1-pAW[i]
seAW[i]=sqrt(pAW[i]*qAW[i]/n1)
LAWD[i]=pAW[i]-(cv*seAW[i])
UAWD[i]=pAW[i]+(cv*seAW[i])

if(LAWD[i]<0) LABB[i]="YES" else LABB[i]="NO"
if(LAWD[i]<0) LAWD[i]=0

if(UAWD[i]>1) UABB[i]="YES" else UABB[i]="NO"
if(UAWD[i]>1) UAWD[i]=1

if(UAWD[i]-LAWD[i]==0)ZWI[i]="YES" else ZWI[i]="NO"
}
return(data.frame(x,LAWD,UAWD,LABB,UABB,ZWI))
}
##############################################################################
#' Adjusted Score method of CI estimation
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - adding factor
#' @details  A score test approach is used after the given data
#' \code{x} and \code{n} are modified as \eqn{x + h} and \eqn{n + (2*h)}
#'  respectively, where \eqn{h > 0} and for all \eqn{x = 0, 1, 2 ..n.}
#' @return A dataframe with
#'  \item{x}{  Number of successes (positive samples)}
#'  \item{LASC }{   Adjusted Score Lower limit}
#'  \item{UASC }{   Adjusted Score Upper Limit}
#'  \item{LABB }{   Adjusted Score Lower Abberation}
#'  \item{UABB }{   Adjusted Score Upper Abberation}
#'  \item{ZWI }{   Zero Width Interval}
#' @family Adjusted methods of CI estimation
#' @seealso \code{\link{prop.test} and \link{binom.test}} for equivalent base Stats R functionality,
#'    \code{\link[binom]{binom.confint}}  provides similar functionality for 11 methods,
#'    \code{\link[PropCIs]{wald2ci}} which provides multiple functions for CI calculation ,
#'    \code{\link[BlakerCI]{binom.blaker.limits}} which calculates Blaker CI which is not covered here and
#'    \code{\link[prevalence]{propCI}} which provides similar functionality.
#' @examples
#' n=5; alp=0.05;h=2
#' ciASC(n,alp,h)
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
ciASC<-function(n,alp,h)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h) >1|| h<0  ) stop("'h' has to be greater than or equal to 0")

####INPUT n
x=0:n
k=n+1
y=x+h
n1=n+(2*h)
####INITIALIZATIONS
pAS=0
qAS=0
seAS=0
LASC=0
UASC=0
LABB=0
UABB=0
ZWI=0

###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
cv1=(cv^2)/(2*n1)
cv2=(cv/(2*n1))^2

#ASCORE (WILSON) METHOD
for(i in 1:k)
{
pAS[i]=y[i]/n1
qAS[i]=1-pAS[i]
seAS[i]=sqrt((pAS[i]*qAS[i]/n1)+cv2)
LASC[i]=(n1/(n1+(cv)^2))*((pAS[i]+cv1)-(cv*seAS[i]))
UASC[i]=(n1/(n1+(cv)^2))*((pAS[i]+cv1)+(cv*seAS[i]))

if(LASC[i]<0) LABB[i]="YES" else LABB[i]="NO"
if(LASC[i]<0) LASC[i]=0

if(UASC[i]>1) UABB[i]="YES" else UABB[i]="NO"
if(UASC[i]>1) UASC[i]=1

if(UASC[i]-LASC[i]==0)ZWI[i]="YES" else ZWI[i]="NO"
}
return(data.frame(x,LASC,UASC,LABB,UABB,ZWI))
}
##############################################################################
##############################################################################
#' Adjusted ArcSine method of CI estimation
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - adding factor
#' @details  Wald-type interval for the arcsine transformation of the parameter
#' \code{p} for the modified data \eqn{x + h} and \eqn{n + (2*h)} , where
#' \eqn{h > 0} and for all \eqn{x = 0, 1, 2 ..n.}
#' @return A dataframe with
#'  \item{x}{  Number of successes (positive samples)}
#'  \item{LAAS }{   Adjusted ArcSine Lower limit}
#'  \item{UAAS }{   Adjusted ArcSine Upper Limit}
#'  \item{LABB }{   Adjusted ArcSine Lower Abberation}
#'  \item{UABB }{   Adjusted ArcSine Upper Abberation}
#'  \item{ZWI }{   Zero Width Interval}
#' @family Adjusted methods of CI estimation
#' @seealso \code{\link{prop.test} and \link{binom.test}} for equivalent base Stats R functionality,
#'    \code{\link[binom]{binom.confint}}  provides similar functionality for 11 methods,
#'    \code{\link[PropCIs]{wald2ci}} which provides multiple functions for CI calculation ,
#'    \code{\link[BlakerCI]{binom.blaker.limits}} which calculates Blaker CI which is not covered here and
#'    \code{\link[prevalence]{propCI}} which provides similar functionality.
#' @examples
#' n=5; alp=0.05;h=2
#' ciAAS(n,alp,h)
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
ciAAS<-function(n,alp,h)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h) >1|| h<0  ) stop("'h' has to be greater than or equal to 0")

####INPUT
x=0:n
k=n+1
y=x+h
m=n+(2*h)
####INITIALIZATIONS
pA=0
qA=0
seA=0
LAAS=0
UAAS=0
LABB=0
UABB=0
ZWI=0

cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#ARC-SINE METHOD
for(i in 1:k)
{
pA[i]=y[i]/m
qA[i]=1-pA[i]
seA[i]=cv/sqrt(4*m)
LAAS[i]=(sin(asin(sqrt(pA[i]))-seA[i]))^2
UAAS[i]=(sin(asin(sqrt(pA[i]))+seA[i]))^2

if(LAAS[i]<0) LABB[i]="YES" else LABB[i]="NO"
if(LAAS[i]<0) LAAS[i]=0

if(UAAS[i]>1) UABB[i]="YES" else UABB[i]="NO"
if(UAAS[i]>1) UAAS[i]=1

if(UAAS[i]-LAAS[i]==0)ZWI[i]="YES" else ZWI[i]="NO"
}
return(data.frame(x,LAAS,UAAS,LABB,UABB,ZWI))
}
##############################################################################
#' Adjusted Likelihood method of CI estimation
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - adding factor
#' @details  Likelihood ratio limits for the data \eqn{x + h} and \eqn{n + (2*h)}
#'  instead of the given code{x} and \code{n}, where \code{h} is a positive integer
#'  \eqn{(1, 2.)} and for all \eqn{x = 0, 1, 2 ..n.}
#' @return A dataframe with
#'  \item{x}{  Number of successes (positive samples)}
#'  \item{LALR }{   Adjusted Likelihood Lower limit}
#'  \item{UALR }{   Adjusted Likelihood Upper Limit}
#'  \item{LABB }{   Adjusted Likelihood Lower Abberation}
#'  \item{UABB }{   Adjusted Likelihood Upper Abberation}
#'  \item{ZWI }{   Zero Width Interval}
#' @family Adjusted methods of CI estimation
#' @seealso \code{\link{prop.test} and \link{binom.test}} for equivalent base Stats R functionality,
#'    \code{\link[binom]{binom.confint}}  provides similar functionality for 11 methods,
#'    \code{\link[PropCIs]{wald2ci}} which provides multiple functions for CI calculation ,
#'    \code{\link[BlakerCI]{binom.blaker.limits}} which calculates Blaker CI which is not covered here and
#'    \code{\link[prevalence]{propCI}} which provides similar functionality.
#' @examples
#' n=5; alp=0.05;h=2
#' ciALR(n,alp,h)
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
ciALR<-function(n,alp,h)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h) >1|| h<0 || !(h%%1 ==0) ) stop("'h' has to be and integer greater than or equal to 0")

####INPUT n
x=0:n
y1=x+h
k=n+1
n1=n+(2*h)
####INITIALIZATIONS
mle=0
cutoff=0
LALR=0
UALR=0
LABB=0
UABB=0
ZWI=0

###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#LIKELIHOOD-RATIO METHOD
for(i in 1:k)
{
likelhd = function(p) dbinom(y1[i],n1,p)
loglik = function(p) dbinom(y1[i],n1,p,log=TRUE)
mle[i]=optimize(likelhd,c(0,1),maximum=TRUE)$maximum
cutoff[i]=loglik(mle[i])-(cv^2/2)
loglik.optim=function(p){abs(cutoff[i]-loglik(p))}
LALR[i]=optimize(loglik.optim, c(0,mle[i]))$minimum
UALR[i]=optimize(loglik.optim, c(mle[i],1))$minimum

if(LALR[i]<0) LABB[i]="YES" else LABB[i]="NO"

if(UALR[i]>1) UABB[i]="YES" else UABB[i]="NO"

if(UALR[i]-LALR[i]==0)ZWI[i]="YES" else ZWI[i]="NO"

}
return(data.frame(x,LALR,UALR,LABB,UABB,ZWI))
}

##############################################################################
#' Adjusted  WALD-T method of CI estimation
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - adding factor
#' @details  Given data \code{x} and \code{n} are modified as \eqn{x + h} and \eqn{n + (2*h)}
#'  respectively, where \eqn{h > 0} then approximate method based on a t_approximation of
#'  the standardized point estimator for all \eqn{x = 0, 1, 2 ..n.}
#' @return A dataframe with
#'  \item{x}{  Number of successes (positive samples)}
#'  \item{LATW }{   Adjusted  WALD-T Lower limit}
#'  \item{UATW }{   Adjusted  WALD-T Upper Limit}
#'  \item{LABB }{   Adjusted  WALD-T Lower Abberation}
#'  \item{UABB }{   Adjusted  WALD-T Upper Abberation}
#'  \item{ZWI }{   Zero Width Interval}
#' @family Adjusted methods of CI estimation
#' @seealso \code{\link{prop.test} and \link{binom.test}} for equivalent base Stats R functionality,
#'    \code{\link[binom]{binom.confint}}  provides similar functionality for 11 methods,
#'    \code{\link[PropCIs]{wald2ci}} which provides multiple functions for CI calculation ,
#'    \code{\link[BlakerCI]{binom.blaker.limits}} which calculates Blaker CI which is not covered here and
#'    \code{\link[prevalence]{propCI}} which provides similar functionality.
#' @examples
#' n=5; alp=0.05;h=2
#' ciATW(n,alp,h)
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
ciATW<-function(n,alp,h)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h) >1|| h<0  ) stop("'h' has to be greater than or equal to 0")

####INPUT n
x=0:n
k=n+1
y=x+h
n1=n+(2*h)
####INITIALIZATIONS
pATW=0
qATW=0
seATW=0
LATW=0
UATW=0
DOF=0
cv=0
LABB=0
UABB=0
ZWI=0
								#Coverage probabilty
#MODIFIED_t-ADJ_WALD METHOD
for(i in 1:k)
{
pATW[i]=y[i]/n1
qATW[i]=1-pATW[i]
f1=function(p,n) p*(1-p)/n
f2=function(p,n) (p*(1-p)/(n^3))+(p+((6*n)-7)*(p^2)+(4*(n-1)*(n-3)*(p^3))-(2*(n-1)*((2*n)-3)*(p^4)))/(n^5)-(2*(p+((2*n)-3)*(p^2)-2*(n-1)*(p^3)))/(n^4)
DOF[i]=2*((f1(pATW[i],n1))^2)/f2(pATW[i],n1)
cv[i]=qt(1-(alp/2), df=DOF[i])
seATW[i]=cv[i]*sqrt(f1(pATW[i],n1))
LATW[i]=pATW[i]-(seATW[i])
UATW[i]=pATW[i]+(seATW[i])

if(LATW[i]<0) LABB[i]="YES" else LABB[i]="NO"
if(LATW[i]<0) LATW[i]=0

if(UATW[i]>1) UABB[i]="YES" else  UABB[i]="NO"
if(UATW[i]>1) UATW[i]=1

if(UATW[i]-LATW[i]==0)ZWI[i]="YES" else ZWI[i]="NO"
}
return(data.frame(x,LATW,UATW,LABB,UABB,ZWI))
}
#####################################################################
##############################################################################
#' Adjusted Logit-Wald method of CI estimation
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - adding factor
#' @details  Wald-type interval for the logit transformation \eqn{log(p/1-p)}
#'  of the parameter \code{p} for the modified data \eqn{x + h} and \eqn{n + (2*h)},
#'  where \eqn{h > 0} and for all \eqn{x = 0, 1, 2 ..n.}
#' @return A dataframe with
#'  \item{x}{  Number of successes (positive samples)}
#'  \item{LALT }{   Adjusted Logit-Wald Lower limit}
#'  \item{UALT }{   Adjusted Logit-Wald Upper Limit}
#'  \item{LABB }{   Adjusted Logit-Wald Lower Abberation}
#'  \item{UABB }{   Adjusted Logit-Wald Upper Abberation}
#'  \item{ZWI }{   Zero Width Interval}
#' @family Adjusted methods of CI estimation
#' @seealso \code{\link{prop.test} and \link{binom.test}} for equivalent base Stats R functionality,
#'    \code{\link[binom]{binom.confint}}  provides similar functionality for 11 methods,
#'    \code{\link[PropCIs]{wald2ci}} which provides multiple functions for CI calculation ,
#'    \code{\link[BlakerCI]{binom.blaker.limits}} which calculates Blaker CI which is not covered here and
#'    \code{\link[prevalence]{propCI}} which provides similar functionality.
#' @examples
#' n=5; alp=0.05;h=2
#' ciALT(n,alp,h)
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
ciALT<-function(n,alp,h)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h) >1|| h<0  ) stop("'h' has to be greater than or equal to 0")

####INPUT n
x=0:n
k=n+1
y=x+h
n1=n+(2*h)

####INITIALIZATIONS
pALT=0
qALT=0
seALT=0
lgit=0
LALT=0
UALT=0
LABB=0
UABB=0
ZWI=0
###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#LOGIT-WALD METHOD
for(i in 1:k)
{
pALT[i]=y[i]/n1
qALT[i]=1-pALT[i]
lgit[i]=log(pALT[i]/qALT[i])
seALT[i]=sqrt(pALT[i]*qALT[i]*n1)
LALT[i]=1/(1+exp(-lgit[i]+(cv/seALT[i])))
UALT[i]=1/(1+exp(-lgit[i]-(cv/seALT[i])))

if(LALT[i]<0) LABB[i]="YES" else LABB[i]="NO"
if(LALT[i]<0) LALT[i]=0

if(UALT[i]>1) UABB[i]="YES" else  UABB[i]="NO"
if(UALT[i]>1) UALT[i]=1

if(UALT[i]-LALT[i]==0)ZWI[i]="YES" else ZWI[i]="NO"
}
return(data.frame(x,LALT,UALT,LABB,UABB,ZWI))
}

#####################################################################
#' CI estimation of 6 adjusted methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)  given adding factor
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - adding factor
#' @details  The Confidence Interval using 6 adjusted methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine) for \code{n} given \code{alp} and \code{h}
#' @return A dataframe with
##'  \item{name }{- Name of the method}
##'  \item{x }{- Number of successes (positive samples)}
##'  \item{LLT }{ - Lower limit}
##'  \item{ULT }{ - Upper Limit}
##'  \item{LABB }{ - Lower Abberation}
##'  \item{UABB }{ - Upper Abberation}
##'  \item{ZWI }{ - Zero Width Interval}
#' @family Adjusted methods of CI estimation
#' @seealso \code{\link{prop.test} and \link{binom.test}} for equivalent base Stats R functionality,
#'    \code{\link[binom]{binom.confint}}  provides similar functionality for 11 methods,
#'    \code{\link[PropCIs]{wald2ci}} which provides multiple functions for CI calculation ,
#'    \code{\link[BlakerCI]{binom.blaker.limits}} which calculates Blaker CI which is not covered here and
#'    \code{\link[prevalence]{propCI}} which provides similar functionality.
#' @examples
#' n=5; alp=0.05;h=2
#' ciAAll(n,alp,h)
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
#7.All methods
ciAAll<-function(n,alp,h)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h) >1|| h<0  || !(h%%1 ==0)) stop("'h' has to be an integer greater than or equal to 0")

  #### Calling functions and creating df
  WaldCI.df    = ciAWD(n,alp,h)
  ArcSineCI.df = ciAAS(n,alp,h)
  LRCI.df      = ciALR(n,alp,round(h,0))		#h must be +ve integer
  ScoreCI.df   = ciASC(n,alp,h)
  WaldLCI.df   = ciALT(n,alp,h)
  AdWaldCI.df  = ciATW(n,alp,h)

  WaldCI.df$method    = as.factor("Adj-Wald")
  ArcSineCI.df$method = as.factor("Adj-ArcSine")
  LRCI.df$method      = as.factor("Adj-Likelihood")
  WaldLCI.df$method    = as.factor("Adj-Logit-Wald")
  ScoreCI.df$method   = as.factor("Adj-Score")
  AdWaldCI.df$method  = as.factor("Adj-Wald-T")

  Generic.1 = data.frame(method = WaldCI.df$method, x=WaldCI.df$x, LowerLimit = WaldCI.df$LAWD, UpperLimit = WaldCI.df$UAWD, LowerAbb = WaldCI.df$LABB, UpperAbb = WaldCI.df$UABB, ZWI = WaldCI.df$ZWI)
  Generic.2 = data.frame(method = ArcSineCI.df$method, x=ArcSineCI.df$x, LowerLimit = ArcSineCI.df$LAAS, UpperLimit = ArcSineCI.df$UAAS, LowerAbb = ArcSineCI.df$LABB, UpperAbb = ArcSineCI.df$UABB, ZWI = ArcSineCI.df$ZWI)
  Generic.3 = data.frame(method = LRCI.df$method, x=LRCI.df$x, LowerLimit = LRCI.df$LALR, UpperLimit = LRCI.df$UALR, LowerAbb = LRCI.df$LABB, UpperAbb = LRCI.df$UABB, ZWI = LRCI.df$ZWI)
  Generic.4 = data.frame(method = ScoreCI.df$method, x=ScoreCI.df$x, LowerLimit = ScoreCI.df$LASC, UpperLimit = ScoreCI.df$UASC, LowerAbb = ScoreCI.df$LABB, UpperAbb = ScoreCI.df$UABB, ZWI = ScoreCI.df$ZWI)
  Generic.5 = data.frame(method = WaldLCI.df$method, x=WaldLCI.df$x, LowerLimit = WaldLCI.df$LALT, UpperLimit = WaldLCI.df$UALT, LowerAbb = WaldLCI.df$LABB, UpperAbb = WaldLCI.df$UABB, ZWI = WaldLCI.df$ZWI)
  Generic.6 = data.frame(method = AdWaldCI.df$method, x=AdWaldCI.df$x, LowerLimit = AdWaldCI.df$LATW, UpperLimit = AdWaldCI.df$UATW, LowerAbb = AdWaldCI.df$LABB, UpperAbb = AdWaldCI.df$UABB, ZWI = AdWaldCI.df$ZWI)

  Final.df= rbind(Generic.1,Generic.2,Generic.3,Generic.4,Generic.5, Generic.6)

  return(Final.df)
}

