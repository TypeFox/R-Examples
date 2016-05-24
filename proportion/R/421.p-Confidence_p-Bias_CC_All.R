#' p-Confidence and p-Bias estimation for continuity corrected Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correction
#' @details  Evaluation of Wald-type interval with
#' continuity correction using p-confidence and p-bias for the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{x1}{  Number of successes (positive samples)}
#'  \item{pconf }{   p-Confidence}
#'  \item{pbias }{   p-Bias}
#' @family p-confidence and p-bias of continuity corrected methods
#' @examples
#' n=5; alp=0.05;c=1/(2*n)
#' pCOpBICWD(n,alp,c)
#' @references
#' [1] 2005 Vos PW and Hudson S.
#' Evaluation Criteria for Discrete Confidence Intervals: Beyond Coverage and Length.
#' The American Statistician: 59; 137 - 142.
#' @export
##### 1.CC-WALD- p-confidence and p-bias for a given n and alpha level
pCOpBICWD<-function(n,alp,c)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
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
pcon=0						#p-confidence
pconC=0
pconf=0
pbia1=0					#p-bias
pbias=0
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
if(LCW[i]<0) LCW[i]=0
if(UCW[i]>1) UCW[i]=1
}
####p-confidence and p-bias
for(i in 2:(k-1))
{
pcon[i-1]=2*(pbinom(i-1, n, LCW[i], lower.tail = FALSE, log.p = FALSE)+dbinom(i-1, n, LCW[i]))
pconC[i-1]=2*pbinom(i-1, n, UCW[i], lower.tail = TRUE, log.p = FALSE)
pconf[i-1]=(1-max(pcon[i-1],pconC[i-1]))*100 		#p-confidence calculation
pbia1[i-1]=max(pcon[i-1],pconC[i-1])-min(pcon[i-1],pconC[i-1])
pbias[i-1]=max(0,pbia1[i-1])*100
}
x1=1:(n-1)
p_C_B=data.frame(x1,pconf,pbias)
return(p_C_B)
}
######################################################################################
#' p-Confidence and p-Bias estimation for continuity corrected Score method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correction
#' @details  Evaluation of continuity corrected score test approach
#' using p-confidence and p-bias for  the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{x1}{  Number of successes (positive samples)}
#'  \item{pconf }{   p-Confidence}
#'  \item{pbias }{   p-Bias}
#' @family p-confidence and p-bias of continuity corrected methods
#' @examples
#' n=5; alp=0.05;c=1/(2*n)
#' pCOpBICSC(n,alp,c)
#' @references
#' [1] 2005 Vos PW and Hudson S.
#' Evaluation Criteria for Discrete Confidence Intervals: Beyond Coverage and Length.
#' The American Statistician: 59; 137 - 142.
#' @export
##### 2.CC-SCORE- p-confidence and p-bias for a given n and alpha level
pCOpBICSC<-function(n,alp,c)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if (c<=0 || c>(1/(2*n)) || length(c)>1) stop("'c' has to be positive and less than or equal to 1/(2*n)")


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
pcon=0						#p-confidence
pconC=0
pconf=0
pbia1=0					#p-bias
pbias=0
###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
cv1=(cv^2)/(2*n)
cv2=cv/(2*n)

#SCORE (WILSON) METHOD
for(i in 1:k)
{
pCS[i]=x[i]/n
qCS[i]=1-pCS[i]
seCS_L[i]=sqrt((cv^2)-(4*n*(c+c^2))+(4*n*pCS[i]*(1-pCS[i]+(2*c))))	#Sq. root term of LL
seCS_U[i]=sqrt((cv^2)+(4*n*(c-c^2))+(4*n*pCS[i]*(1-pCS[i]-(2*c))))	#Sq. root term of LL
LCS[i]=(n/(n+(cv)^2))*((pCS[i]-c+cv1)-(cv2*seCS_L[i]))
UCS[i]=(n/(n+(cv)^2))*((pCS[i]+c+cv1)+(cv2*seCS_U[i]))
if(LCS[i]<0) LCS[i]=0
if(UCS[i]>1) UCS[i]=1
}
####p-confidence and p-bias
for(i in 2:(k-1))
{
pcon[i-1]=2*(pbinom(i-1, n, LCS[i], lower.tail = FALSE, log.p = FALSE)+dbinom(i-1, n, LCS[i]))
pconC[i-1]=2*pbinom(i-1, n, UCS[i], lower.tail = TRUE, log.p = FALSE)
pconf[i-1]=1-max(pcon[i-1],pconC[i-1]) 		#p-confidence calculation
pbia1[i-1]=max(pcon[i-1],pconC[i-1])-min(pcon[i-1],pconC[i-1])
pbias[i-1]=max(0,pbia1[i-1])
}
x1=1:(n-1)
p_C_B=data.frame(x1,pconf,pbias)
return(p_C_B)
}
#################################################################################
#' p-Confidence and p-Bias estimation for continuity corrected ArcSine method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correction
#' @details  Evaluation of continuity corrected Wald-type interval for the arcsine
#' transformation of the parameter \code{p} using p-confidence and p-bias for the  \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{x1}{  Number of successes (positive samples)}
#'  \item{pconf }{   p-Confidence}
#'  \item{pbias }{   p-Bias}
#' @family p-confidence and p-bias of continuity corrected methods
#' @examples
#' n=5; alp=0.05;c=1/(2*n)
#' pCOpBICAS(n,alp,c)
#' @references
#' [1] 2005 Vos PW and Hudson S.
#' Evaluation Criteria for Discrete Confidence Intervals: Beyond Coverage and Length.
#' The American Statistician: 59; 137 - 142.
#' @export
##### 3.CC-ARC SINE - p-confidence and p-bias for a given n and alpha level
pCOpBICAS<-function(n,alp,c)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
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
pCA[i]=x[i]/n
qCA[i]=1-pCA[i]
seCA[i]=cv/sqrt(4*n)
LCA[i]=(sin(asin(sqrt(pCA[i]))-seCA[i]-c))^2
UCA[i]=(sin(asin(sqrt(pCA[i]))+seCA[i]+c))^2
if(LCA[i]<0) LCA[i]=0
if(UCA[i]>1) UCA[i]=1
}

####p-confidence and p-bias
for(i in 2:(k-1))
{
  pcon[i-1]=2*(pbinom(i-1, n, LCA[i], lower.tail = FALSE, log.p = FALSE)+dbinom(i-1, n, LCA[i]))
  pconC[i-1]=2*pbinom(i-1, n, UCA[i], lower.tail = TRUE, log.p = FALSE)
  pconf[i-1]=1-max(pcon[i-1],pconC[i-1]) 		#p-confidence calculation
  pbia1[i-1]=max(pcon[i-1],pconC[i-1])-min(pcon[i-1],pconC[i-1])
  pbias[i-1]=as.numeric(max(0,pbia1[i-1]))
}

x1=1:(n-1)
p_C_B=data.frame(x1,pconf,pbias)
return(p_C_B)
}

#########################################################################################
#' p-Confidence and p-Bias estimation for continuity corrected Logit Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correction
#' @details  Evaluation of continuity corrected Wald-type interval based on the logit
#' transformation of \code{p} using p-confidence and p-bias for the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{x1}{  Number of successes (positive samples)}
#'  \item{pconf }{   p-Confidence}
#'  \item{pbias }{   p-Bias}
#' @family p-confidence and p-bias of continuity corrected methods
#' @examples
#' n=5; alp=0.05;c=1/(2*n)
#' pCOpBICLT(n,alp,c)
#' @references
#' [1] 2005 Vos PW and Hudson S.
#' Evaluation Criteria for Discrete Confidence Intervals: Beyond Coverage and Length.
#' The American Statistician: 59; 137 - 142.
#' @export
##### 4.CC - LOGIT WALD - p-confidence and p-bias for a given n and alpha level
pCOpBICLT<-function(n,alp,c)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
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
pcon=0						#p-confidence
pconC=0
pconf=0
pbia1=0					#p-bias
pbias=0
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
if(LCLT[j+1]<0) LCLT[j+1]=0
if(UCLT[j+1]>1) UCLT[j+1]=1
}

####p-confidence and p-bias
for(i in 2:(k-1))
{
pcon[i-1]=2*(pbinom(i-1, n, LCLT[i], lower.tail = FALSE, log.p = FALSE)+dbinom(i-1, n, LCLT[i]))
pconC[i-1]=2*pbinom(i-1, n, UCLT[i], lower.tail = TRUE, log.p = FALSE)
pconf[i-1]=1-max(pcon[i-1],pconC[i-1]) 		#p-confidence calculation
pbia1[i-1]=max(pcon[i-1],pconC[i-1])-min(pcon[i-1],pconC[i-1])
pbias[i-1]=as.numeric(max(0,pbia1[i-1]))
}
x1=1:(n-1)
p_C_B=data.frame(x1,pconf,pbias)
return(p_C_B)
}

####################################################################################
#' p-Confidence and p-Bias estimation for continuity corrected Wald-T method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correction
#' @details  Evaluation of continuity corrected Wald method based on a t_approximation
#' of the standardized point estimator using p-confidence and p-bias for the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{x1}{  Number of successes (positive samples)}
#'  \item{pconf }{   p-Confidence}
#'  \item{pbias }{   p-Bias}
#' @family p-confidence and p-bias of continuity corrected methods
#' @examples
#' n=5; alp=0.05;c=1/(2*n)
#' pCOpBICTW(n,alp,c)
#' @references
#' [1] 2005 Vos PW and Hudson S.
#' Evaluation Criteria for Discrete Confidence Intervals: Beyond Coverage and Length.
#' The American Statistician: 59; 137 - 142.
#' @export
##### 5.T-WALD_CC: p-confidence and p-bias for a given n and alpha level
pCOpBICTW<-function(n,alp,c)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
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
pcon=0						#p-confidence
pconC=0
pconf=0
pbia1=0					#p-bias
pbias=0
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
if(LCTW[i]<0) LCTW[i]=0
if(UCTW[i]>1) UCTW[i]=1
}
####p-confidence and p-bias
for(i in 2:(k-1))
{
pcon[i-1]=2*(pbinom(i-1, n, LCTW[i], lower.tail = FALSE, log.p = FALSE)+dbinom(i-1, n, LCTW[i]))
pconC[i-1]=2*pbinom(i-1, n, UCTW[i], lower.tail = TRUE, log.p = FALSE)
pconf[i-1]=1-max(pcon[i-1],pconC[i-1]) 		#p-confidence calculation
pbia1[i-1]=max(pcon[i-1],pconC[i-1])-min(pcon[i-1],pconC[i-1])
pbias[i-1]=as.numeric(max(0,pbia1[i-1]))
}
x1=1:(n-1)
p_C_B=data.frame(x1,pconf,pbias)
return(p_C_B)
}

#####################################################################################################################################
#' Performs p-Confidence and p-Bias estimation for 5 continuity corrected methods (Wald, Wald-T,  Score, Logit-Wald, ArcSine)
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correction
#' @details  Evaluation of p-Confidence and p-Bias estimation of 5 continuity corrected methods (Wald, Wald-T,  Score, Logit-Wald, ArcSine)
#' @return A dataframe with
#'  \item{x1}{  Number of successes (positive samples)}
#'  \item{pconf }{   p-Confidence}
#'  \item{pbias }{   p-Bias}
#'  \item{method}{   Method name}
#' @family p-confidence and p-bias of continuity corrected methods
#' @examples
#' n=5; alp=0.05;c=1/(2*n)
#' pCOpBICAll(n,alp,c)
#' @references
#' [1] 2005 Vos PW and Hudson S.
#' Evaluation Criteria for Discrete Confidence Intervals: Beyond Coverage and Length.
#' The American Statistician: 59; 137 - 142.
#' @export
#7.All methods
pCOpBICAll<-function(n,alp,c)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if (c<=0 || c>(1/(2*n)) || length(c)>1) stop("'c' has to be positive and less than or equal to 1/(2*n)")

  #### Calling functions and creating df
  WaldpCB.df    = pCOpBICWD(n,alp,c)
  ArcSinepCB.df = pCOpBICAS(n,alp,c)
  ScorepCB.df   = pCOpBICSC(n,alp,c)
  WaldLpCB.df   = pCOpBICLT(n,alp,c)
  AdWaldpCB.df  = pCOpBICTW(n,alp,c)

  WaldpCB.df$method    = as.factor("CC-Wald")
  ArcSinepCB.df$method = as.factor("CC-ArcSine")
  WaldLpCB.df$method   = as.factor("CC-Logit-Wald")
  ScorepCB.df$method   = as.factor("CC-Score")
  AdWaldpCB.df$method  = as.factor("CC-Wald-T")

  Final.df= rbind(WaldpCB.df, ArcSinepCB.df, ScorepCB.df,WaldLpCB.df,AdWaldpCB.df)

  return(Final.df)
}

