#' Performs expected length and sum of length of Adjusted Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  Evaluation of adjusted Wald-type interval using sum of length of the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{sumLen}{  The sum of the expected length}
#'  \item{explMean}{  The mean of the expected length}
#'  \item{explSD}{  The Standard Deviation of the expected length}
#'  \item{explMax}{  The max of the expected length}
#'  \item{explLL}{  The Lower limit of the expected length calculated using mean - SD}
#'  \item{explUL}{  The Upper limit of the expected length calculated using mean + SD}
#' @family Expected length  of adjusted methods
#' @examples
#' n= 10; alp=0.05; h=2;a=1;b=1;
#' lengthAWD(n,alp,h,a,b)
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
##### 1.ADJUSTED WALD sum of length for a given n and alpha level
lengthAWD<-function(n,alp,h,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(a)>1 || h<0  ) stop("'h' has to be greater than or equal to 0")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")


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
s=2000
LEAW=0 								#LENGTH OF INTERVAL

ewiAW=matrix(0,k,s)						#sum of length quantity in sum
ewAW=0									#sum of length
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
LEAW[i]=UAW[i]-LAW[i]
}
#sumLEAW=sum(LEAW)
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
  for(i in 1:k)
  {
    ewiAW[i,j]=LEAW[i]*dbinom(i-1, n,hp[j])
  }
  ewAW[j]=sum(ewiAW[,j])						#Expected Length
}

sumLen=sum(LEAW)
explMean=mean(ewAW)
explSD=sd(ewAW)
explMax=max(ewAW)
explLL=explMean-(explSD)
explUL=explMean+(explSD)
df.Summary=data.frame(sumLen,explMean,explSD,explMax,explLL,explUL)
return(df.Summary)
}
###############################################################################################################
#' Performs expected length and sum of length of Adjusted Score method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  Evaluation of adjusted score test approach using sum of length of the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{sumLen}{  The sum of the expected length}
#'  \item{explMean}{  The mean of the expected length}
#'  \item{explSD}{  The Standard Deviation of the expected length}
#'  \item{explMax}{  The max of the expected length}
#'  \item{explLL}{  The Lower limit of the expected length calculated using mean - SD}
#'  \item{explUL}{  The Upper limit of the expected length calculated using mean + SD}
#' @family Expected length  of adjusted methods
#' @examples
#' n= 10; alp=0.05; h=2;a=1;b=1;
#' lengthASC(n,alp,h,a,b)
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
##### 2.ADJUSTED SCORE - sum of length for a given n and alpha level
lengthASC<-function(n,alp,h,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(a)>1 || h<0  ) stop("'h' has to be greater than or equal to 0")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")

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
s=2000
LEAS=0 								#LENGTH OF INTERVAL

ewiAS=matrix(0,k,s)						#sum of length quantity in sum
ewAS=0									#sum of length
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
LEAS[i]=UAS[i]-LAS[i]
}
#sumLEAS=sum(LEAS)
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
  for(i in 1:k)
  {
    ewiAS[i,j]=LEAS[i]*dbinom(i-1, n,hp[j])
  }
  ewAS[j]=sum(ewiAS[,j])						#Expected Length
}

sumLen=sum(LEAS)
explMean=mean(ewAS)
explSD=sd(ewAS)
explMax=max(ewAS)
explLL=explMean-(explSD)
explUL=explMean+(explSD)
df.Summary=data.frame(sumLen,explMean,explSD,explMax,explLL,explUL)
return(df.Summary)
}

###############################################################################################################
#' Expected length and sum of length of Adjusted ArcSine method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  Evaluation of adjusted Wald-type interval for the arcsine transformation of the parameter p using sum of length of the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{sumLen}{  The sum of the expected length}
#'  \item{explMean}{  The mean of the expected length}
#'  \item{explSD}{  The Standard Deviation of the expected length}
#'  \item{explMax}{  The max of the expected length}
#'  \item{explLL}{  The Lower limit of the expected length calculated using mean - SD}
#'  \item{explUL}{  The Upper limit of the expected length calculated using mean + SD}
#' @family Expected length  of adjusted methods
#' @examples
#' n= 10; alp=0.05; h=2;a=1;b=1;
#' lengthAAS(n,alp,h,a,b)
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
##### 3. ADJUSTED ARC SINE - sum of length for a given n and alpha level
lengthAAS<-function(n,alp,h,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(a)>1 || h<0  ) stop("'h' has to be greater than or equal to 0")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")

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
s=2000
LEAA=0 								#LENGTH OF INTERVAL

ewiAA=matrix(0,k,s)						#sum of length quantity in sum
ewAA=0									#sum of length
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
LEAA[i]=UAA[i]-LAA[i]
}
#sumLEAA=sum(LEAA)
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
  for(i in 1:k)
  {
    ewiAA[i,j]=LEAA[i]*dbinom(i-1, n,hp[j])
  }
  ewAA[j]=sum(ewiAA[,j])						#Expected Length
}

sumLen=sum(LEAA)
explMean=mean(ewAA)
explSD=sd(ewAA)
explMax=max(ewAA)
explLL=explMean-(explSD)
explUL=explMean+(explSD)
df.Summary=data.frame(sumLen,explMean,explSD,explMax,explLL,explUL)
return(df.Summary)
}

###############################################################################################################
#' Performs expected length and sum of length of Adjusted Logit Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  Evaluation of adjusted Wald-type interval based on the logit transformation of p using sum of length of the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{sumLen}{  The sum of the expected length}
#'  \item{explMean}{  The mean of the expected length}
#'  \item{explSD}{  The Standard Deviation of the expected length}
#'  \item{explMax}{  The max of the expected length}
#'  \item{explLL}{  The Lower limit of the expected length calculated using mean - SD}
#'  \item{explUL}{  The Upper limit of the expected length calculated using mean + SD}
#' @family Expected length  of adjusted methods
#' @examples
#' n= 10; alp=0.05; h=2;a=1;b=1;
#' lengthALT(n,alp,h,a,b)
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
##### 4.ADJUSTED LOGIT-WALD - sum of length for a given n and alpha level
lengthALT<-function(n,alp,h,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(a)>1 || h<0  ) stop("'h' has to be greater than or equal to 0")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")

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
s=2000
LEALT=0 								#LENGTH OF INTERVAL

ewiALT=matrix(0,k,s)						#sum of length quantity in sum
ewALT=0									#sum of length
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
LEALT[i]=UALT[i]-LALT[i]
}
#sumLEALT=sum(LEALT)
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
  for(i in 1:k)
  {
    ewiALT[i,j]=LEALT[i]*dbinom(i-1, n,hp[j])
  }
  ewALT[j]=sum(ewiALT[,j])						#Expected Length
}
sumLen=sum(LEALT)
explMean=mean(ewALT)
explSD=sd(ewALT)
explMax=max(ewALT)
explLL=explMean-(explSD)
explUL=explMean+(explSD)
df.Summary=data.frame(sumLen,explMean,explSD,explMax,explLL,explUL)
return(df.Summary)
}
###############################################################################################################
#' Performs expected length and sum of length of Adjusted Wald-T method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  Evaluation of approximate and adjusted method based on a t_approximation of the standardized point estimator using sum of length of the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{sumLen}{  The sum of the expected length}
#'  \item{explMean}{  The mean of the expected length}
#'  \item{explSD}{  The Standard Deviation of the expected length}
#'  \item{explMax}{  The max of the expected length}
#'  \item{explLL}{  The Lower limit of the expected length calculated using mean - SD}
#'  \item{explUL}{  The Upper limit of the expected length calculated using mean + SD}
#' @family Expected length  of adjusted methods
#' @examples
#' n= 10; alp=0.05; h=2;a=1;b=1;
#' lengthATW(n,alp,h,a,b)
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
##### 5.ADJUSTED t-WALD - sum of length for a given n and alpha level
lengthATW<-function(n,alp,h,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(a)>1 || h<0  ) stop("'h' has to be greater than or equal to 0")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")

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
s=2000
LEATW=0 								#LENGTH OF INTERVAL

ewiATW=matrix(0,k,s)						#sum of length quantity in sum
ewATW=0									#sum of length
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
if(LATW[i]<0) LATW[i]=0
if(UATW[i]>1) UATW[i]=1
LEATW[i]=UATW[i]-LATW[i]
}
#sumLEATW=sum(LEATW)
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
  for(i in 1:k)
  {
    ewiATW[i,j]=LEATW[i]*dbinom(i-1, n,hp[j])
  }
  ewATW[j]=sum(ewiATW[,j])						#Expected Length
}

sumLen=sum(LEATW)
explMean=mean(ewATW)
explSD=sd(ewATW)
explMax=max(ewATW)
explLL=explMean-(explSD)
explUL=explMean+(explSD)
df.Summary=data.frame(sumLen,explMean,explSD,explMax,explLL,explUL)
return(df.Summary)
}

###############################################################################################################
#' Performs expected length and sum of length of Adjusted Likelihood method
#' Performs expected length and sum of length of Adjusted Likelihood method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  Evaluation of adjusted Likelihood ratio limits using sum of length of the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{sumLen}{  The sum of the expected length}
#'  \item{explMean}{  The mean of the expected length}
#'  \item{explSD}{  The Standard Deviation of the expected length}
#'  \item{explMax}{  The max of the expected length}
#'  \item{explLL}{  The Lower limit of the expected length calculated using mean - SD}
#'  \item{explUL}{  The Upper limit of the expected length calculated using mean + SD}
#' @family Expected length  of adjusted methods
#' @examples
#' n= 10; alp=0.05; h=2;a=1;b=1;
#' lengthALR(n,alp,h,a,b)
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
#####6.ADJUSTED LIKELIHOOD RATIO - sum of length for a given n and alpha level
lengthALR<-function(n,alp,h,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h) >1|| h<0  || !(h%%1 ==0)) stop("'h' has to be an integer greater than or equal to 0")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")

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
s=2000
LEAL=0 								#LENGTH OF INTERVAL

ewiAL=matrix(0,k,s)						#sum of length quantity in sum
ewAL=0									#sum of length
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
LEAL[i]=UAL[i]-LAL[i]
}
#sumLEAL=sum(LEAL)
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
  for(i in 1:k)
  {
    ewiAL[i,j]=LEAL[i]*dbinom(i-1, n,hp[j])
  }
  ewAL[j]=sum(ewiAL[,j])						#Expected Length
}

sumLen=sum(LEAL)
explMean=mean(ewAL)
explSD=sd(ewAL)
explMax=max(ewAL)
explLL=explMean-(explSD)
explUL=explMean+(explSD)
df.Summary=data.frame(sumLen,explMean,explSD,explMax,explLL,explUL)
return(df.Summary)
}

########################################################################################
#' Expected Length summary calculation using 6 adjusted methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  The sum of length of 6 adjusted methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)  for \code{n} given \code{alp}, \code{a}, \code{b}
#' @return A dataframe with
#'  \item{sumLen}{  The sum of the expected length}
#'  \item{explMean}{  The mean of the expected length}
#'  \item{explSD}{  The Standard Deviation of the expected length}
#'  \item{explMax}{  The max of the expected length}
#'  \item{explLL}{  The Lower limit of the expected length calculated using mean - SD}
#'  \item{explUL}{  The Upper limit of the expected length calculated using mean + SD}
#'  \item{method}{Name of the method}
#' @family Expected length  of adjusted methods
#' @examples
#' n= 10; alp=0.05; h=2;a=1;b=1;
#' lengthAAll(n,alp,h,a,b)
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
lengthAAll<-function(n,alp,h,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h) >1|| h<0  || !(h%%1 ==0)) stop("'h' has to be an integer greater than or equal to 0")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")

  #### Calling functions and creating df

  df1    = lengthAWD(n,alp,h,a,b)
  df2    = lengthASC(n,alp,h,a,b)
  df3    = lengthAAS(n,alp,h,a,b)
  df4    = lengthALT(n,alp,h,a,b)
  df5    = lengthATW(n,alp,h,a,b)
  df6    = lengthALR(n,alp,h,a,b)

  df1$method = "Adj-Wald"
  df2$method = "Adj-Score"
  df3$method = "Adj-ArcSine"
  df4$method = "Adj-Logit-Wald"
  df5$method = "Adj-Wald-T"
  df6$method = "Adj-Likelihood"

  Final.df= rbind(df1,df2,df3,df4,df5,df6)

  return(Final.df)
}

########################################################################################
# Performs expected length calculation using 6 adjusted methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)
##### 10.All methods - Expected length
explAAll<-function(n,alp,h,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(a)>1 || h<0  ) stop("'h' has to be greater than or equal to 0")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")

  #### Calling functions and creating df
  df.1    = gexplAWD(n,alp,h,a,b)
  df.2    = gexplASC(n,alp,h,a,b)
  df.3    = gexplAAS(n,alp,h,a,b)
  df.4    = gexplALT(n,alp,h,a,b)
  df.5    = gexplATW(n,alp,h,a,b)
  df.6    = gexplALR(n,alp,h,a,b)

  df.new=  rbind(df.1,df.2,df.3,df.4,df.5,df.6)
  return(df.new)
}

