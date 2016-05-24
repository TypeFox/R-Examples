#' p-Confidence and p-Bias estimation for adjusted Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @details  Evaluation of adjusted Wald-type
#' interval using p-confidence and p-bias for the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{x1}{  Number of successes (positive samples)}
#'  \item{pconf }{   p-Confidence}
#'  \item{pbias }{   p-Bias}
#' @family p-confidence and p-bias of adjusted methods
#' @examples
#' n=5; alp=0.05;h=2
#' pCOpBIAWD(n,alp,h)
#' @references
#' [1] 2005 Vos PW and Hudson S.
#' Evaluation Criteria for Discrete Confidence Intervals: Beyond Coverage and Length.
#' The American Statistician: 59; 137 - 142.
#' @export
##### 1.ADJUSTED WALD- p-confidence and p-bias for a given n and alpha level
pCOpBIAWD<-function(n,alp,h)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if ((class(h) != "integer") & (class(h) != "numeric") || h<0  ) stop("'h' has to be greater than or equal to 0")


####INPUT n
x=0:n
k=n+1
y=x+h
m=n+(2*h)
####INITIALIZATIONS
pAW=0
qAW=0
seAW=0
LAW=0
UAW=0
pcon=0					#p-confidence
pconC=0
pconf=0
pbia1=0					#p-bias
pbias=0
###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#WALD METHOD
for(i in 1:k)
{
pAW[i]=y[i]/m
qAW[i]=1-pAW[i]
seAW[i]=sqrt(pAW[i]*qAW[i]/m)
LAW[i]=pAW[i]-(cv*seAW[i])
UAW[i]=pAW[i]+(cv*seAW[i])
if(LAW[i]<0) LAW[i]=0
if(UAW[i]>1) UAW[i]=1
}
####p-confidence and p-bias
for(i in 2:(k-1))
{
pcon[i-1]=2*(pbinom(i-1, n, LAW[i], lower.tail = FALSE, log.p = FALSE)+dbinom(i-1, n, LAW[i]))
pconC[i-1]=2*pbinom(i-1, n, UAW[i], lower.tail = TRUE, log.p = FALSE)
pconf[i-1]=(1-max(pcon[i-1],pconC[i-1]))*100 		#p-confidence calculation
pbia1[i-1]=max(pcon[i-1],pconC[i-1])-min(pcon[i-1],pconC[i-1])
pbias[i-1]=max(0,pbia1[i-1])*100
}
x1=1:(n-1)
p_C_B=data.frame(x1,pconf,pbias)
return(p_C_B)
}

#####################################################################################################################################
#' p-Confidence and p-Bias estimation for adjusted Score method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @details  Evaluation of adjusted score test approach
#' using p-confidence and p-bias for the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{x1}{  Number of successes (positive samples)}
#'  \item{pconf }{   p-Confidence}
#'  \item{pbias }{   p-Bias}
#' @family p-confidence and p-bias of adjusted methods
#' @examples
#' n=5; alp=0.05;h=2
#' pCOpBIASC(n,alp,h)
#' @references
#' [1] 2005 Vos PW and Hudson S.
#' Evaluation Criteria for Discrete Confidence Intervals: Beyond Coverage and Length.
#' The American Statistician: 59; 137 - 142.
#' @export
##### 2.ADJUSTED SCORE- p-confidence and p-bias for a given n and alpha level
pCOpBIASC<-function(n,alp,h)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if ((class(h) != "integer") & (class(h) != "numeric") || h<0  ) stop("'h' has to be greater than or equal to 0")

####INPUT n
x=0:n
k=n+1
y=x+h
n1=n+(2*h)
####INITIALIZATIONS
pAS=0
qAS=0
seAS=0
LAS=0
UAS=0
pcon=0						#p-confidence
pconC=0
pconf=0
pbia1=0					#p-bias
pbias=0
###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
cv1=(cv^2)/(2*n1)
cv2=(cv/(2*n1))^2

#SCORE (WILSON) METHOD
for(i in 1:k)
{
pAS[i]=y[i]/n1
qAS[i]=1-pAS[i]
seAS[i]=sqrt((pAS[i]*qAS[i]/n1)+cv2)
LAS[i]=(n1/(n1+(cv)^2))*((pAS[i]+cv1)-(cv*seAS[i]))
UAS[i]=(n1/(n1+(cv)^2))*((pAS[i]+cv1)+(cv*seAS[i]))
if(LAS[i]<0) LAS[i]=0
if(UAS[i]>1) UAS[i]=1
}
####p-confidence and p-bias
for(i in 2:(k-1))
{
pcon[i-1]=2*(pbinom(i-1, n, LAS[i], lower.tail = FALSE, log.p = FALSE)+dbinom(i-1, n, LAS[i]))
pconC[i-1]=2*pbinom(i-1, n, UAS[i], lower.tail = TRUE, log.p = FALSE)
pconf[i-1]=1-max(pcon[i-1],pconC[i-1]) 		#p-confidence calculation
pbia1[i-1]=max(pcon[i-1],pconC[i-1])-min(pcon[i-1],pconC[i-1])
pbias[i-1]=as.numeric(max(0,pbia1[i-1]))
}
x1=1:(n-1)
p_C_B=data.frame(x1,pconf,pbias)
return(p_C_B)
}
#####################################################################################################################################
#' p-Confidence and p-Bias estimation for adjusted ArcSine method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @details  Evaluation of adjusted Wald-type interval for the arcsine
#' transformation of the parameter \code{p} using p-confidence and p-bias for the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{x1}{  Number of successes (positive samples)}
#'  \item{pconf }{   p-Confidence}
#'  \item{pbias }{   p-Bias}
#' @family p-confidence and p-bias of adjusted methods
#' @examples
#' n=5; alp=0.05;h=2
#' pCOpBIAAS(n,alp,h)
#' @references
#' [1] 2005 Vos PW and Hudson S.
#' Evaluation Criteria for Discrete Confidence Intervals: Beyond Coverage and Length.
#' The American Statistician: 59; 137 - 142.
#' @export
##### 3.ADJUSTED ARC SINE - p-confidence and p-bias for a given n and alpha level
pCOpBIAAS<-function(n,alp,h)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if ((class(h) != "integer") & (class(h) != "numeric") || h<0  ) stop("'h' has to be greater than or equal to 0")

####INPUT n
x=0:n
k=n+1
y=x+h
m=n+(2*h)
####INITIALIZATIONS
pAA=0
qAA=0
seAA=0
LAA=0
UAA=0
pcon=0						#p-confidence
pconC=0
pconf=0
pbia1=0					#p-bias
pbias=0
###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#ARC-SINE METHOD
for(i in 1:k)
{
pAA[i]=y[i]/m
qAA[i]=1-pAA[i]
seAA[i]=cv/sqrt(4*m)
LAA[i]=(sin(asin(sqrt(pAA[i]))-seAA[i]))^2
UAA[i]=(sin(asin(sqrt(pAA[i]))+seAA[i]))^2
if(LAA[i]<0) LAA[i]=0
if(UAA[i]>1) UAA[i]=1
}

####p-confidence and p-bias
for(i in 2:(k-1))
{
pcon[i-1]=2*(pbinom(i-1, n, LAA[i], lower.tail = FALSE, log.p = FALSE)+dbinom(i-1, n, LAA[i]))
pconC[i-1]=2*pbinom(i-1, n, UAA[i], lower.tail = TRUE, log.p = FALSE)
pconf[i-1]=1-max(pcon[i-1],pconC[i-1]) 		#p-confidence calculation
pbia1[i-1]=max(pcon[i-1],pconC[i-1])-min(pcon[i-1],pconC[i-1])
pbias[i-1]=as.numeric(max(0,pbia1[i-1]))
}

x1=1:(n-1)
p_C_B=data.frame(x1,pconf,pbias)
return(p_C_B)
}
#####################################################################################################################################
#' p-Confidence and p-Bias estimation for adjusted Logit-Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @details  Evaluation of adjusted Wald-type interval based on the logit transformation of
#' \code{p} using p-confidence and p-bias for the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{x1}{  Number of successes (positive samples)}
#'  \item{pconf }{   p-Confidence}
#'  \item{pbias }{   p-Bias}
#' @family p-confidence and p-bias of adjusted methods
#' @examples
#' n=5; alp=0.05;h=2
#' pCOpBIALT(n,alp,h)
#' @references
#' [1] 2005 Vos PW and Hudson S.
#' Evaluation Criteria for Discrete Confidence Intervals: Beyond Coverage and Length.
#' The American Statistician: 59; 137 - 142.
#' @export
##### 4.ADJUSTED LOGIT WALD - p-confidence and p-bias for a given n and alpha level
pCOpBIALT<-function(n,alp,h)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if ((class(h) != "integer") & (class(h) != "numeric") || h<0  ) stop("'h' has to be greater than or equal to 0")

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
pcon=0						#p-confidence
pconC=0
pconf=0
pbia1=0					#p-bias
pbias=0
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
if(LALT[i]<0) LALT[i]=0
if(UALT[i]>1) UALT[i]=1
}

####p-confidence and p-bias
for(i in 2:(k-1))
{
pcon[i-1]=2*(pbinom(i-1, n, LALT[i], lower.tail = FALSE, log.p = FALSE)+dbinom(i-1, n, LALT[i]))
pconC[i-1]=2*pbinom(i-1, n, UALT[i], lower.tail = TRUE, log.p = FALSE)
pconf[i-1]=1-max(pcon[i-1],pconC[i-1]) 		#p-confidence calculation
pbia1[i-1]=max(pcon[i-1],pconC[i-1])-min(pcon[i-1],pconC[i-1])
pbias[i-1]=as.numeric(max(0,pbia1[i-1]))
}
x1=1:(n-1)
p_C_B=data.frame(x1,pconf,pbias)
return(p_C_B)
}

#####################################################################################################################################
#' p-Confidence and p-Bias estimation for adjusted  Wald-T method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @details  Evaluation of approximate Wald method based on a t_approximation
#' of the standardized point estimator using p-confidence and p-bias for the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{x1}{  Number of successes (positive samples)}
#'  \item{pconf }{   p-Confidence}
#'  \item{pbias }{   p-Bias}
#' @family p-confidence and p-bias of adjusted methods
#' @examples
#' n=5; alp=0.05;h=2
#' pCOpBIATW(n,alp,h)
#' @references
#' [1] 2005 Vos PW and Hudson S.
#' Evaluation Criteria for Discrete Confidence Intervals: Beyond Coverage and Length.
#' The American Statistician: 59; 137 - 142.
#' @export
##### 5.ADJUSTED T-WALD: p-confidence and p-bias for a given n and alpha level
pCOpBIATW<-function(n,alp,h)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if ((class(h) != "integer") & (class(h) != "numeric") || h<0  ) stop("'h' has to be greater than or equal to 0")

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
pcon=0						#p-confidence
pconC=0
pconf=0
pbia1=0					#p-bias
pbias=0
#MODIFIED_t-WALD METHOD
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
if(LATW[i]<0) LATW[i]=0
if(UATW[i]>1) UATW[i]=1
}
####p-confidence and p-bias
for(i in 2:(k-1))
{
pcon[i-1]=2*(pbinom(i-1, n, LATW[i], lower.tail = FALSE, log.p = FALSE)+dbinom(i-1, n, LATW[i]))
pconC[i-1]=2*pbinom(i-1, n, UATW[i], lower.tail = TRUE, log.p = FALSE)
pconf[i-1]=1-max(pcon[i-1],pconC[i-1]) 		#p-confidence calculation
pbia1[i-1]=max(pcon[i-1],pconC[i-1])-min(pcon[i-1],pconC[i-1])
pbias[i-1]=as.numeric(max(0,pbia1[i-1]))
}
x1=1:(n-1)
p_C_B=data.frame(x1,pconf,pbias)
return(p_C_B)
}

###########################################################################################
#' p-Confidence and p-Bias estimation for adjusted Likelihood method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @details  Evaluation of adjusted Likelihood ratio limits
#' using p-confidence and p-bias for the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{x1}{  Number of successes (positive samples)}
#'  \item{pconf }{   p-Confidence}
#'  \item{pbias }{   p-Bias}
#' @family p-confidence and p-bias of adjusted methods
#' @examples
#' n=5; alp=0.05;h=2
#' pCOpBIALR(n,alp,h)
#' @references
#' [1] 2005 Vos PW and Hudson S.
#' Evaluation Criteria for Discrete Confidence Intervals: Beyond Coverage and Length.
#' The American Statistician: 59; 137 - 142.
#' @export
##### 6.ADJUSTED LIKELIHOOD RATIO - p-confidence and p-bias for a given n and alpha level
pCOpBIALR<-function(n,alp,h)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h) >1|| h<0  || !(h%%1 ==0)) stop("'h' has to be an integer greater than or equal to 0")

####INPUT n
y=0:n
y1=y+h
k=n+1
n1=n+(2*h)
####INITIALIZATIONS
mle=0
cutoff=0
LAL=0
UAL=0

pcon=0						#p-confidence
pconC=0
pconf=0
pbia1=0					#p-bias
pbias=0


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
LAL[i]=optimize(loglik.optim, c(0,mle[i]))$minimum
UAL[i]=optimize(loglik.optim, c(mle[i],1))$minimum
}
####p-confidence and p-bias
for(i in 2:(k-1))
{
pcon[i-1]=2*(pbinom(i-1, n, LAL[i], lower.tail = FALSE, log.p = FALSE)+dbinom(i-1, n, LAL[i]))
pconC[i-1]=2*pbinom(i-1, n, UAL[i], lower.tail = TRUE, log.p = FALSE)
pconf[i-1]=1-max(pcon[i-1],pconC[i-1]) 		#p-confidence calculation
pbia1[i-1]=max(pcon[i-1],pconC[i-1])-min(pcon[i-1],pconC[i-1])
pbias[i-1]=as.numeric(max(0,pbia1[i-1]))
}
x1=1:(n-1)
p_C_B=data.frame(x1,pconf,pbias)
return(p_C_B)
}
#####################################################################################################################################
#'  Performs p-Confidence and p-Bias estimation of 6 adjusted methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine) given adding factor
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @details  Evaluation of p-Confidence and p-Bias estimation of 6 adjusted methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)
#' @return A dataframe with
#'  \item{x1}{  Number of successes (positive samples)}
#'  \item{pconf }{   p-Confidence}
#'  \item{pbias }{   p-Bias}
#'  \item{method}{Method name}
#' @family p-confidence and p-bias of adjusted methods
#' @examples
#' n=5; alp=0.05;h=2
#' pCOpBIAAll(n,alp,h)
#' @references
#' [1] 2005 Vos PW and Hudson S.
#' Evaluation Criteria for Discrete Confidence Intervals: Beyond Coverage and Length.
#' The American Statistician: 59; 137 - 142.
#' @export
#7.All methods
pCOpBIAAll<-function(n,alp,h)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h) >1|| h<0  || !(h%%1 ==0)) stop("'h' has to be an integer greater than or equal to 0")

  #### Calling functions and creating df
  WaldpCB.df    = pCOpBIAWD(n,alp,h)
  ArcSinepCB.df = pCOpBIAAS(n,alp,h)
  LRpCB.df      = pCOpBIALR(n,alp,h)
  ScorepCB.df   = pCOpBIASC(n,alp,h)
  WaldLpCB.df   = pCOpBIALT(n,alp,h)
  AdWaldpCB.df  = pCOpBIATW(n,alp,h)

  WaldpCB.df$method    = as.factor("Adj-Wald")
  ArcSinepCB.df$method = as.factor("Adj-ArcSine")
  LRpCB.df$method      = as.factor("Adj-Likelihood")
  WaldLpCB.df$method   = as.factor("Adj-Logit-Wald")
  ScorepCB.df$method   = as.factor("Adj-Score")
  AdWaldpCB.df$method  = as.factor("Adj-Wald-T")

  Final.df= rbind(WaldpCB.df, ArcSinepCB.df, LRpCB.df,ScorepCB.df,WaldLpCB.df,AdWaldpCB.df)

  return(Final.df)
}
