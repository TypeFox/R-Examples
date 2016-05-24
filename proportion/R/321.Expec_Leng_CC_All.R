#' Expected Length summary of continuity corrected Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correction
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  Evaluation of Wald-type interval with continuity correction using sum of length of the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{sumLen}{  The sum of the expected length}
#'  \item{explMean}{  The mean of the expected length}
#'  \item{explSD}{  The Standard Deviation of the expected length}
#'  \item{explMax}{  The max of the expected length}
#'  \item{explLL}{  The Lower limit of the expected length calculated using mean - SD}
#'  \item{explUL}{  The Upper limit of the expected length calculated using mean + SD}
#' @family Expected length  of continuity corrected methods
#' @examples
#' n= 10; alp=0.05; c=1/(2*n);a=1;b=1;
#' lengthCWD(n,alp,c,a,b)
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
##### 1.CC-WALD sum of length for a given n and alpha level
lengthCWD<-function(n,alp,c,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(c) != "integer") & (class(c) != "numeric") || length(c) >1 || c<0 ) stop("'c' has to be positive")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")

####INPUT n
x=0:n
k=n+1
####INITIALIZATIONS
pCW=0
qCW=0
seCW=0
LCW=0
UCW=0
s=5000
LECW=0 								#LENGTH OF INTERVAL

ewiCW=matrix(0,k,s)						#sum of length quantity in sum
ewCW=0									#sum of length
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
LECW[i]=UCW[i]-LCW[i]
}
#sumLECW=sum(LECW)
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
  for(i in 1:k)
  {
    ewiCW[i,j]=LECW[i]*dbinom(i-1, n,hp[j])
  }
  ewCW[j]=sum(ewiCW[,j])						#Expected Length
}

sumLen=sum(LECW)
explMean=mean(ewCW)
explSD=sd(ewCW)
explMax=max(ewCW)
explLL=explMean-(explSD)
explUL=explMean+(explSD)
df.length=data.frame(sumLen,explMean,explSD,explMax,explLL,explUL)
return(df.length)
}
###############################################################################################################
#' Expected Length summary of continuity corrected Score method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correction
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  Evaluation of continuity corrected score test approach using sum of length of the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{sumLen}{  The sum of the expected length}
#'  \item{explMean}{  The mean of the expected length}
#'  \item{explSD}{  The Standard Deviation of the expected length}
#'  \item{explMax}{  The max of the expected length}
#'  \item{explLL}{  The Lower limit of the expected length calculated using mean - SD}
#'  \item{explUL}{  The Upper limit of the expected length calculated using mean + SD}
#' @family Expected length  of continuity corrected methods
#' @examples
#' n= 10; alp=0.05; c=1/(2*n);a=1;b=1;
#' lengthCSC(n,alp,c,a,b)
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
##### 2.CC SCORE - sum of length for a given n and alpha level
lengthCSC<-function(n,alp,c,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(c) != "integer") & (class(c) != "numeric") || length(c) >1 || c<0 ) stop("'c' has to be positive")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")

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
s=5000
LECS=0 								#LENGTH OF INTERVAL

ewiCS=matrix(0,k,s)						#sum of length quantity in sum
ewCS=0									#sum of length
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
if(LCS[i]<0) LCS[i]=0
if(UCS[i]>1) UCS[i]=1
LECS[i]=UCS[i]-LCS[i]
}
#sumLECS=sum(LECS)
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
  for(i in 1:k)
  {
    ewiCS[i,j]=LECS[i]*dbinom(i-1, n,hp[j])
  }
  ewCS[j]=sum(ewiCS[,j])						#Expected Length
}

sumLen=sum(LECS)
explMean=mean(ewCS)
explSD=sd(ewCS)
explMax=max(ewCS)
explLL=explMean-(explSD)
explUL=explMean+(explSD)
df.length=data.frame(sumLen,explMean,explSD,explMax,explLL,explUL)
return(df.length)
}
###############################################################################################################
#' Expected Length summary of continuity corrected ArcSine method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correction
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  Evaluation of continuity corrected Wald-type interval for the arcsine transformation of the parameter p using sum of length of the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{sumLen}{  The sum of the expected length}
#'  \item{explMean}{  The mean of the expected length}
#'  \item{explSD}{  The Standard Deviation of the expected length}
#'  \item{explMax}{  The max of the expected length}
#'  \item{explLL}{  The Lower limit of the expected length calculated using mean - SD}
#'  \item{explUL}{  The Upper limit of the expected length calculated using mean + SD}
#' @family Expected length  of continuity corrected methods
#' @examples
#' n= 10; alp=0.05; c=1/(2*n);a=1;b=1;
#' lengthCAS(n,alp,c,a,b)
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
##### 3.CC ARC SINE - sum of length for a given n and alpha level
lengthCAS<-function(n,alp,c,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(c) != "integer") & (class(c) != "numeric") || length(c) >1 || c<0 ) stop("'c' has to be positive")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")

####INPUT n
x=0:n
k=n+1
####INITIALIZATIONS
pCA=0
qCA=0
seCA=0
LCA=0
UCA=0
s=5000
LECA=0 								#LENGTH OF INTERVAL

ewiCA=matrix(0,k,s)						#sum of length quantity in sum
ewCA=0									#sum of length
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
LECA[i]=UCA[i]-LCA[i]
}
#sumLECA=sum(LECA)
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
  for(i in 1:k)
  {
    ewiCA[i,j]=LECA[i]*dbinom(i-1, n,hp[j])
  }
  ewCA[j]=sum(ewiCA[,j])						#Expected Length
}

sumLen=sum(LECA)
explMean=mean(ewCA)
explSD=sd(ewCA)
explMax=max(ewCA)
explLL=explMean-(explSD)
explUL=explMean+(explSD)
df.length=data.frame(sumLen,explMean,explSD,explMax,explLL,explUL)
return(df.length)
}

###############################################################################################################
#' Expected Length summary of continuity corrected Logit Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correction
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  Evaluation of continuity corrected Wald-type interval based on the logit transformation of p using sum of length of the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{sumLen}{  The sum of the expected length}
#'  \item{explMean}{  The mean of the expected length}
#'  \item{explSD}{  The Standard Deviation of the expected length}
#'  \item{explMax}{  The max of the expected length}
#'  \item{explLL}{  The Lower limit of the expected length calculated using mean - SD}
#'  \item{explUL}{  The Upper limit of the expected length calculated using mean + SD}
#' @family Expected length  of continuity corrected methods
#' @examples
#' n= 10; alp=0.05; c=1/(2*n);a=1;b=1;
#' lengthCLT(n,alp,c,a,b)
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
##### 4.CC LOGIT-WALD - sum of length for a given n and alpha level
lengthCLT<-function(n,alp,c,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(c) != "integer") & (class(c) != "numeric") || length(c) >1 || c<0 ) stop("'c' has to be positive")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")

####INPUT n
x=0:n
k=n+1
								#sum of length
#INITIALIZATIONS
pCLT=0
qCLT=0
seCLT=0
lgit=0
LCLT=0
UCLT=0
s=5000
LECLT=0 								#LENGTH OF INTERVAL

ewiCLT=matrix(0,k,s)						#sum of length quantity in sum
ewCLT=0									#sum of length
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
if(LCLT[i]<0) LCLT[i]=0
if(UCLT[i]>1) UCLT[i]=1
LECLT[i]=UCLT[i]-LCLT[i]
}
#sumLECLT=sum(LECLT)
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
  for(i in 1:k)
  {
    ewiCLT[i,j]=LECLT[i]*dbinom(i-1, n,hp[j])
  }
  ewCLT[j]=sum(ewiCLT[,j])						#Expected Length
}

sumLen=sum(LECLT)
explMean=mean(ewCLT)
explSD=sd(ewCLT)
explMax=max(ewCLT)
explLL=explMean-(explSD)
explUL=explMean+(explSD)
df.length=data.frame(sumLen,explMean,explSD,explMax,explLL,explUL)
return(df.length)
}

###############################################################################################################
#' Expected Length summary of continuity corrected Wald-T method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correction
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  Evaluation of approximate and continuity corrected method based on
#' a t_approximation of the standardized point estimator using sum of length
#' of the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{sumLen}{  The sum of the expected length}
#'  \item{explMean}{  The mean of the expected length}
#'  \item{explSD}{  The Standard Deviation of the expected length}
#'  \item{explMax}{  The max of the expected length}
#'  \item{explLL}{  The Lower limit of the expected length calculated using mean - SD}
#'  \item{explUL}{  The Upper limit of the expected length calculated using mean + SD}
#' @family Expected length  of continuity corrected methods
#' @examples
#' n= 10; alp=0.05; c=1/(2*n);a=1;b=1;
#' lengthCTW(n,alp,c,a,b)
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
##### 5.CC t-WALD_CC - sum of length for a given n and alpha level
lengthCTW<-function(n,alp,c,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(c) != "integer") & (class(c) != "numeric") || length(c) >1 || c<0 ) stop("'c' has to be positive")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")

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
s=5000
LECTW=0 								#LENGTH OF INTERVAL

ewiCTW=matrix(0,k,s)						#sum of length quantity in sum
ewCTW=0									#sum of length
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
LECTW[i]=UCTW[i]-LCTW[i]
}
#sumLECTW=sum(LECTW)
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
  for(i in 1:k)
  {
    ewiCTW[i,j]=LECTW[i]*dbinom(i-1, n,hp[j])
  }
  ewCTW[j]=sum(ewiCTW[,j])						#Expected Length
}

sumLen=sum(LECTW)
explMean=mean(ewCTW)
explSD=sd(ewCTW)
explMax=max(ewCTW)
explLL=explMean-(explSD)
explUL=explMean+(explSD)
df.length=data.frame(sumLen,explMean,explSD,explMax,explLL,explUL)
return(df.length)
}
########################################################################################
#' Expected Length summary calculation using 5 continuity corrected methods
#' (Wald, Wald-T, Score, Logit-Wald, ArcSine)
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correction
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  The sum of length of 5 continuity corrected methods (Wald, Wald-T, Score, Logit-Wald, ArcSine) of \code{n} given \code{alp}, \code{a}, \code{b}
#' @return A dataframe with
#'  \item{sumLen}{  The sum of the expected length}
#'  \item{explMean}{  The mean of the expected length}
#'  \item{explSD}{  The Standard Deviation of the expected length}
#'  \item{explMax}{  The max of the expected length}
#'  \item{explLL}{  The Lower limit of the expected length calculated using mean - SD}
#'  \item{explUL}{  The Upper limit of the expected length calculated using mean + SD}
#' @family Expected length  of continuity corrected methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; c=1/(2*n);a=1;b=1;
#' lengthCAll(n,alp,c,a,b)
#' }
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
##### 9. sum of length for a given n and alpha level for all methods
lengthCAll<-function(n,alp,c,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if (c<=0 || c>(1/(2*n)) || length(c)>1) stop("'c' has to be positive and less than or equal to 1/(2*n)")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")

  #### Calling functions and creating df

  df1    = lengthCWD(n,alp,c,a,b)
  df2    = lengthCSC(n,alp,c,a,b)
  df3    = lengthCAS(n,alp,c,a,b)
  df4    = lengthCLT(n,alp,c,a,b)
  df5    = lengthCTW(n,alp,c,a,b)

  df1$method = "CC-Wald"
  df2$method = "CC-Score"
  df3$method = "CC-ArcSine"
  df4$method = "CC-Logit-Wald"
  df5$method = "CC-Wald-T"

  Final.df= rbind(df1,df2,df3,df4,df5)

  return(Final.df)
}

########################################################################################
# Expected length calculation using 5 continuity corrected methods (Wald, Wald-T, Score, Logit-Wald, ArcSine)
##### 10.. Expected length for a given n and alpha level for all methods
explCAll<-function(n,alp,c,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if (c<=0 || c>(1/(2*n)) || length(c)>1) stop("'c' has to be positive and less than or equal to 1/(2*n)")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")

  #### Calling functions and creating df
  df.1    = gexplCWD(n,alp,c,a,b)
  df.2    = gexplCSC(n,alp,c,a,b)
  df.3    = gexplCAS(n,alp,c,a,b)
  df.4    = gexplCLT(n,alp,c,a,b)
  df.5    = gexplCTW(n,alp,c,a,b)

  df.new=  rbind(df.1,df.2,df.3,df.4,df.5)
  return(df.new)
}
