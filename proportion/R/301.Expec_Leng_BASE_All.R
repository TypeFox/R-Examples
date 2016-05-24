#' Expected length and sum of length of  Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  Evaluation of Wald-type intervals using sum of length of the \eqn{n + 1}
#'  intervals
#' @return A dataframe with
#'  \item{sumLen}{  The sum of the expected length}
#'  \item{explMean}{  The mean of the expected length}
#'  \item{explSD}{  The Standard Deviation of the expected length}
#'  \item{explMax}{  The max of the expected length}
#'  \item{explLL}{  The Lower limit of the expected length calculated using mean - SD}
#'  \item{explUL}{  The Upper limit of the expected length calculated using mean + SD}
#' @family Expected length  of base methods
#' @examples
#' n=5; alp=0.05;a=1;b=1
#' lengthWD(n,alp,a,b)
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
##### 1.WALD sum of length for a given n and alpha level
lengthWD<-function(n,alp,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")

####INPUT n
x=0:n
k=n+1
####INITIALIZATIONS
pW=0
qW=0
seW=0
LW=0
UW=0
s=5000
LEW=0 								#LENGTH OF INTERVAL

ewiW=matrix(0,k,s)						#sum of length quantity in sum
ewW=0									#sum of length
###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#WALD METHOD
for(i in 1:k)
{
pW[i]=x[i]/n
qW[i]=1-(x[i]/n)
seW[i]=sqrt(pW[i]*qW[i]/n)
LW[i]=pW[i]-(cv*seW[i])
UW[i]=pW[i]+(cv*seW[i])
if(LW[i]<0) LW[i]=0
if(UW[i]>1) UW[i]=1
LEW[i]=UW[i]-LW[i]
}
#sumLEW=sum(LEW)
####sum of length
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
ewiW[i,j]=LEW[i]*dbinom(i-1, n,hp[j])
}
ewW[j]=sum(ewiW[,j])						#sum of length
}
#ELW=data.frame(hp,ewW)
sumLen=sum(LEW)
explMean=mean(ewW)
explSD=sd(ewW)
explMax=max(ewW)
explLL=explMean-(explSD)
explUL=explMean+(explSD)
df.Summary=data.frame(sumLen,explMean,explSD,explMax,explLL,explUL)
return(df.Summary)
}
###############################################################################################################
#' Expected length and sum of length of  Score method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  Evaluation of score test approach using sum of length of the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{sumLen}{  The sum of the expected length}
#'  \item{explMean}{  The mean of the expected length}
#'  \item{explSD}{  The Standard Deviation of the expected length}
#'  \item{explMax}{  The max of the expected length}
#'  \item{explLL}{  The Lower limit of the expected length calculated using mean - SD}
#'  \item{explUL}{  The Upper limit of the expected length calculated using mean + SD}
#' @family Expected length  of base methods
#' @examples
#' n=5; alp=0.05;a=1;b=1
#' lengthSC(n,alp,a,b)
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
##### 2.SCORE - sum of length for a given n and alpha level
lengthSC<-function(n,alp,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")

####INPUT n
x=0:n
k=n+1
####INITIALIZATIONS
pS=0
qS=0
seS=0
LS=0
US=0
s=5000
LES=0 								#LENGTH OF INTERVAL
ewiS=matrix(0,k,s)						#sum of length quantity in sum
ewS=0									#sum of length
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
LS[i]=(n/(n+(cv)^2))*((pS[i]+cv1)-(cv*seS[i]))
US[i]=(n/(n+(cv)^2))*((pS[i]+cv1)+(cv*seS[i]))
if(LS[i]<0) LS[i]=0
if(US[i]>1) US[i]=1
LES[i]=US[i]-LS[i]
}
#sumLES=sum(LES)

####sum of length
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
ewiS[i,j]=LES[i]*dbinom(i-1, n,hp[j])
}
ewS[j]=sum(ewiS[,j])						#sum of length
}
#ELS=data.frame(hp,ewS)
sumLen=sum(LES)
explMean=mean(ewS)
explSD=sd(ewS)
explMax=max(ewS)
explLL=explMean-(explSD)
explUL=explMean+(explSD)
df.Summary=data.frame(sumLen,explMean,explSD,explMax,explLL,explUL)
return(df.Summary)
}
###############################################################################################################
#' Expected length and sum of length of  ArcSine method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  Evaluation of Wald-type interval for the arcsine transformation of the parameter
#' \code{p} using sum of length of the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{sumLen}{  The sum of the expected length}
#'  \item{explMean}{  The mean of the expected length}
#'  \item{explSD}{  The Standard Deviation of the expected length}
#'  \item{explMax}{  The max of the expected length}
#'  \item{explLL}{  The Lower limit of the expected length calculated using mean - SD}
#'  \item{explUL}{  The Upper limit of the expected length calculated using mean + SD}
#' @family Expected length  of base methods
#' @examples
#' n=5; alp=0.05;a=1;b=1
#' lengthAS(n,alp,a,b)
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
##### 3. ARC SINE - sum of length for a given n and alpha level
lengthAS<-function(n,alp,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")

####INPUT n
x=0:n
k=n+1
####INITIALIZATIONS
pA=0
qA=0
seA=0
LA=0
UA=0
s=5000
LEA=0 								#LENGTH OF INTERVAL

ewiA=matrix(0,k,s)						#sum of length quantity in sum
ewA=0									#sum of length
###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#WALD METHOD
#ARC-SINE METHOD
for(i in 1:k)
{
pA[i]=x[i]/n
qA[i]=1-pA[i]
seA[i]=cv/sqrt(4*n)
LA[i]=(sin(asin(sqrt(pA[i]))-seA[i]))^2
UA[i]=(sin(asin(sqrt(pA[i]))+seA[i]))^2
if(LA[i]<0) LA[i]=0
if(UA[i]>1) UA[i]=1
LEA[i]=UA[i]-LA[i]
}
#sumLEA=sum(LEA)
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
  for(i in 1:k)
  {
    ewiA[i,j]=LEA[i]*dbinom(i-1, n,hp[j])
  }
  ewA[j]=sum(ewiA[,j])						#Expected Length
}

sumLen=sum(LEA)
explMean=mean(ewA)
explSD=sd(ewA)
explMax=max(ewA)
explLL=explMean-(explSD)
explUL=explMean+(explSD)
df.Summary=data.frame(sumLen,explMean,explSD,explMax,explLL,explUL)
return(df.Summary)
}

###############################################################################################################
#' Expected length and sum of length of  Logit Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  Evaluation of Wald-type interval based on the logit
#' transformation of \code{p} using sum of length of the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{sumLen}{  The sum of the expected length}
#'  \item{explMean}{  The mean of the expected length}
#'  \item{explSD}{  The Standard Deviation of the expected length}
#'  \item{explMax}{  The max of the expected length}
#'  \item{explLL}{  The Lower limit of the expected length calculated using mean - SD}
#'  \item{explUL}{  The Upper limit of the expected length calculated using mean + SD}
#' @family Expected length  of base methods
#' @examples
#' n=5; alp=0.05;a=1;b=1
#' lengthLT(n,alp,a,b)
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
##### 4.LOGIT-WALD - sum of length for a given n and alpha level
lengthLT<-function(n,alp,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")

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
s=5000
LELT=0 								#LENGTH OF INTERVAL

ewiLT=matrix(0,k,s)						#sum of length quantity in sum
ewLT=0									#sum of length
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
if(LLT[j+1]<0) LLT[j+1]=0
if(ULT[j+1]>1) ULT[j+1]=1
}
for(i in 1:k)
{
LELT[i]=ULT[i]-LLT[i]
}
#sumLET=sum(LELT)
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
  for(i in 1:k)
  {
    ewiLT[i,j]=LELT[i]*dbinom(i-1, n,hp[j])
  }
  ewLT[j]=sum(ewiLT[,j])						#Expected Length
}

sumLen=sum(LELT)
explMean=mean(ewLT)
explSD=sd(ewLT)
explMax=max(ewLT)
explLL=explMean-(explSD)
explUL=explMean+(explSD)
df.Summary=data.frame(sumLen,explMean,explSD,explMax,explLL,explUL)
return(df.Summary)
}

###############################################################################################################
#' Expected length and sum of length of  Wald-T method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  Evaluation of approximate method based on a t_approximation of the
#' standardized point estimator using sum of length of the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{sumLen}{  The sum of the expected length}
#'  \item{explMean}{  The mean of the expected length}
#'  \item{explSD}{  The Standard Deviation of the expected length}
#'  \item{explMax}{  The max of the expected length}
#'  \item{explLL}{  The Lower limit of the expected length calculated using mean - SD}
#'  \item{explUL}{  The Upper limit of the expected length calculated using mean + SD}
#' @family Expected length  of base methods
#' @examples
#' n=5; alp=0.05;a=1;b=1
#' lengthTW(n,alp,a,b)
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
##### 5.t-WALD - sum of length for a given n and alpha level
lengthTW<-function(n,alp,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")

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
s=5000
LETW=0 								#LENGTH OF INTERVAL

ewiTW=matrix(0,k,s)						#sum of length quantity in sum
ewTW=0									#sum of length
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
if(LTW[i]<0) LTW[i]=0
if(UTW[i]>1) UTW[i]=1
LETW[i]=UTW[i]-LTW[i]
}
#sumLETW=sum(LETW)
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
  for(i in 1:k)
  {
    ewiTW[i,j]=LETW[i]*dbinom(i-1, n,hp[j])
  }
  ewTW[j]=sum(ewiTW[,j])						#Expected Length
}

sumLen=sum(LETW)
explMean=mean(ewTW)
explSD=sd(ewTW)
explMax=max(ewTW)
explLL=explMean-(explSD)
explUL=explMean+(explSD)
df.Summary=data.frame(sumLen,explMean,explSD,explMax,explLL,explUL)
return(df.Summary)
}

###############################################################################################################
#' Expected length and sum of length of  likelihood Ratio method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  Evaluation of Likelihood ratio limits using sum of length of the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{sumLen}{  The sum of the expected length}
#'  \item{explMean}{  The mean of the expected length}
#'  \item{explSD}{  The Standard Deviation of the expected length}
#'  \item{explMax}{  The max of the expected length}
#'  \item{explLL}{  The Lower limit of the expected length calculated using mean - SD}
#'  \item{explUL}{  The Upper limit of the expected length calculated using mean + SD}
#' @family Expected length  of base methods
#' @examples
#' n=5; alp=0.05;a=1;b=1
#' lengthLR(n,alp,a,b)
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
#####6.LIKELIHOOD RATIO - sum of length for a given n and alpha level
lengthLR<-function(n,alp,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")

####INPUT n
y=0:n
k=n+1
####INITIALIZATIONS
mle=0
cutoff=0
LL=0
UL=0
s=5000
LEL=0 								#LENGTH OF INTERVAL

ewiL=matrix(0,k,s)						#sum of length quantity in sum
ewL=0										#Simulation run to generate hypothetical p
###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#LIKELIHOOD-RATIO METHOD
for(i in 1:k)
{
likelhd = function(p) dbinom(y[i],n,p)
loglik = function(p) dbinom(y[i],n,p,log=TRUE)
mle[i]=optimize(likelhd,c(0,1),maximum=TRUE)$maximum
cutoff[i]=loglik(mle[i])-(cv^2/2)
loglik.optim=function(p){abs(cutoff[i]-loglik(p))}
LL[i]=optimize(loglik.optim, c(0,mle[i]))$minimum
UL[i]=optimize(loglik.optim, c(mle[i],1))$minimum
LEL[i]=UL[i]-LL[i]
}
#sumLEL=sum(LEL)
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
  for(i in 1:k)
  {
    ewiL[i,j]=LEL[i]*dbinom(i-1, n,hp[j])
  }
  ewL[j]=sum(ewiL[,j])						#Expected Length
}

sumLen=sum(LEL)
explMean=mean(ewL)
explSD=sd(ewL)
explMax=max(ewL)
explLL=explMean-(explSD)
explUL=explMean+(explSD)
df.Summary=data.frame(sumLen,explMean,explSD,explMax,explLL,explUL)
return(df.Summary)
}
###################################################################################################
#' Expected length and sum of length of  Exact method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param e - Exact method indicator  in [0, 1] {1: Clopper Pearson, 0.5: Mid P}
#' The input can also be a range of values between 0 and 1.
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  Evaluation of Confidence interval for \code{p} based on inverting equal-tailed
#' binomial tests with null hypothesis \eqn{H0: p = p0} using sum of length of the \eqn{n + 1} intervals.
#' @return A dataframe with
#'  \item{sumLen}{  The sum of the expected length}
#'  \item{explMean}{  The mean of the expected length}
#'  \item{explSD}{  The Standard Deviation of the expected length}
#'  \item{explMax}{  The max of the expected length}
#'  \item{explLL}{  The Lower limit of the expected length calculated using mean - SD}
#'  \item{explUL}{  The Upper limit of the expected length calculated using mean + SD}
#' @family Expected length  of base methods
#' @examples
#' \dontrun{
#' n=5; alp=0.05;e=0.5;a=1;b=1
#' lengthEX(n,alp,e,a,b)
#' n=5; alp=0.05;e=1;a=1;b=1 #Clopper-Pearson
#' lengthEX(n,alp,e,a,b)
#' n=5; alp=0.05;e=c(0.1,0.5,0.95,1);a=1;b=1 #Range including Mid-p and Clopper-Pearson
#' lengthEX(n,alp,e,a,b)
#' }
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
##### 7.EXACT EMTHOD sum of length for a given n and alpha level
lengthEX<-function(n,alp,e,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(e)) stop("'e' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(e) != "integer") & (class(e) != "numeric") || any(e>1) || any(e<0)) stop("'e' has to be between 0 and 1")
  if (length(e)>10) stop("'e' can have only 10 intervals")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")

  nvar=length(e)

  res <- data.frame()

  for(i in 1:nvar)
  {
    lu=glengthEX301(n,alp,e[i],a,b)
    res <- rbind(res,lu)
  }
  return(res)
}

glengthEX301<-function(n,alp,e,a,b)
  {
####INPUT n
  x=0:n
  k=n+1
####INITIALIZATIONS
LEX=0
UEX=0
s=5000
LEEX=0 								#LENGTH OF INTERVAL
ewiEX=matrix(0,k,s)						#sum of length quantity in sum
ewEX=0									#sum of length

#EXACT METHOD
LEX[1]=0
UEX[1]= 1-((alp/(2*e))^(1/n))
LEX[k]=(alp/(2*e))^(1/n)
UEX[k]=1
LEEX[1]=1-((alp/(2*e))^(1/n))
LEEX[k]=1-((alp/(2*e))^(1/n))

for(i in 1:(k-2))
{
LEX[i+1]=exlim301l(x[i+1],n,alp,e)
UEX[i+1]=exlim301u(x[i+1],n,alp,e)
LEEX[i+1]=UEX[i+1]-LEX[i+1]
}
#sumLEEX=sum(LEEX)
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
  for(i in 1:k)
  {
    ewiEX[i,j]=LEEX[i]*dbinom(i-1, n,hp[j])
  }
  ewEX[j]=sum(ewiEX[,j])						#Expected Length
}

sumLen=sum(LEEX)
explMean=mean(ewEX)
explSD=sd(ewEX)
explMax=max(ewEX)
explLL=explMean-(explSD)
explUL=explMean+(explSD)
df.Summary=data.frame(sumLen,explMean,explSD,explMax,explLL,explUL,e)
return(df.Summary)
}
#####TO FIND LOWER LIMITS
exlim301l=function(x,n,alp,e)
{
  z=x-1
  y=0:z
  f1=function(p) (1-e)*dbinom(x,n,p)+sum(dbinom(y,n,p))-(1-(alp/2))
  LEX= uniroot(f1,c(0,1))$root
  return(LEX)
}
#####TO FIND UPPER LIMITS
exlim301u=function(x,n,alp,e)
{
  z=x-1
  y=0:z
  f2  = function(p) e*dbinom(x,n,p)+sum(dbinom(y,n,p))-(alp/2)
  UEX = uniroot(f2,c(0,1))$root
  return(UEX)
}

###############################################################################################################
#' Expected length and sum of length of  Bayesian method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param a1 - Beta Prior Parameters for Bayesian estimation
#' @param a2 - Beta Prior Parameters for Bayesian estimation
#' @details  Evaluation of Bayesian Highest Probability Density (HPD) and two
#' tailed intervals using sum of length of the \eqn{n + 1}
#' intervals for the Beta - Binomial conjugate prior model for the probability of success \code{p}
#' @return A dataframe with
#'  \item{sumLen}{  The sum of the expected length}
#'  \item{explMean}{  The mean of the expected length}
#'  \item{explSD}{  The Standard Deviation of the expected length}
#'  \item{explMax}{  The max of the expected length}
#'  \item{explLL}{  The Lower limit of the expected length calculated using mean - SD}
#'  \item{explUL}{  The Upper limit of the expected length calculated using mean + SD}
#'  \item{method}{  The method used - Quantile and HPD}
#' @family Expected length  of base methods
#' @examples
#' n=5; alp=0.05;a=1;b=1;a1=1;a2=1
#' lengthBA(n,alp,a,b,a1,a2)
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
##### 8.BAYESIAN sum of length for a given n and alpha level
lengthBA<-function(n,alp,a,b,a1,a2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(a1)) stop("'a1' is missing")
  if (missing(a2)) stop("'a2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if ((class(a1) != "integer") & (class(a1) != "numeric") || length(a1)>1 || a1<0  ) stop("'a1' has to be greater than or equal to 0")
  if ((class(a2) != "integer") & (class(a2) != "numeric") || length(a2)>1 || a2<0  ) stop("'a2' has to be greater than or equal to 0")


####INPUT n
x=0:n
k=n+1
####INITIALIZATIONS
LBAQ=0
UBAQ=0
LBAH=0
UBAH=0
s=5000
LEBAQ=0 								#LENGTH OF INTERVAL
LEBAH=0
ewiBAQ=matrix(0,k,s)						#sum of length quantity in sum
ewBAQ=0
ewiBAH=matrix(0,k,s)						#sum of length quantity in sum
ewBAH=0									#sum of length

#library(TeachingDemos)				#To get HPDs
for(i in 1:k)
{
#Quantile Based Intervals
LBAQ[i]=qbeta(alp/2,x[i]+a1,n-x[i]+a2)
UBAQ[i]=qbeta(1-(alp/2),x[i]+a1,n-x[i]+a2)

LBAH[i]=TeachingDemos::hpd(qbeta,shape1=x[i]+a1,shape2=n-x[i]+a2,conf=1-alp)[1]
UBAH[i]=TeachingDemos::hpd(qbeta,shape1=x[i]+a1,shape2=n-x[i]+a2,conf=1-alp)[2]

LEBAQ[i]=UBAQ[i]-LBAQ[i]
LEBAH[i]=UBAH[i]-LBAH[i]
}
# sumLEBAQ=sum(LEBAQ)
# sumLEBAH=sum(LEBAH)
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
  for(i in 1:k)
  {
    ewiBAQ[i,j]=LEBAQ[i]*dbinom(i-1, n,hp[j])
    ewiBAH[i,j]=LEBAH[i]*dbinom(i-1, n,hp[j])

  }
  ewBAQ[j]=sum(ewiBAQ[,j])
  ewBAH[j]=sum(ewiBAH[,j])						#Expected Length
}
sumLenBAQ=sum(LEBAQ)
explMeanBAQ=mean(ewBAQ)
explSDBAQ=sd(ewBAQ)
explMaxBAQ=max(ewBAQ)
explLLBAQ=explMeanBAQ-(explSDBAQ)
explULBAQ=explMeanBAQ+(explSDBAQ)
df.SummaryBAQ=data.frame(sumLen=sumLenBAQ,explMean=explMeanBAQ,
                          explSD=explSDBAQ,explMax=explMaxBAQ,
                          explLL=explLLBAQ,explUL=explULBAQ,method="Quantile")

sumLenBAH=sum(LEBAH)
explMeanBAH=mean(ewBAH)
explSDBAH=sd(ewBAH)
explMaxBAH=max(ewBAH)
explLLBAH=explMeanBAH-(explSDBAH)
explULBAH=explMeanBAH+(explSDBAH)
df.SummaryBAH=data.frame(sumLen=sumLenBAH,explMean=explMeanBAH,
                          explSD=explSDBAH,explMax=explMaxBAH,
                          explLL=explLLBAH,explUL=explULBAH,method="HPD")
df.Summary=rbind(df.SummaryBAQ,df.SummaryBAH)
return(df.Summary)
}
########################################################################################
#' Expected Length summary calculation using 6 base methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  The sum of length for 6 base methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine) for \code{n} given alpha \code{alp} and  Beta parameters \code{a}, \code{b}
#' @return A dataframe with
#'  \item{sumLen}{  The sum of the expected length}
#'  \item{explMean}{  The mean of the expected length}
#'  \item{explSD}{  The Standard Deviation of the expected length}
#'  \item{explMax}{  The max of the expected length}
#'  \item{explLL}{  The Lower limit of the expected length calculated using mean - SD}
#'  \item{explUL}{  The Upper limit of the expected length calculated using mean + SD}
#'  \item{method}{  The name of the method}
#' @family Expected length  of base methods
#' @examples
#' \dontrun{
#' n=5; alp=0.05;a=1;b=1
#' lengthAll(n,alp,a,b)
#' }
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
##### 9. sum of length for a given n and alpha level for all methods
lengthAll<-function(n,alp,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")


  #### Calling functions and creating df

  df1    = lengthWD(n,alp,a,b)
  df2    = lengthSC(n,alp,a,b)
  df3    = lengthAS(n,alp,a,b)
  df4    = lengthLT(n,alp,a,b)
  df5    = lengthTW(n,alp,a,b)
  df6    = lengthLR(n,alp,a,b)

  df1$method = "Wald"
  df2$method = "Score"
  df3$method = "ArcSine"
  df4$method = "Logit-Wald"
  df5$method = "Wald-T"
  df6$method = "Likelihood"

  Final.df= rbind(df1,df2,df3,df4,df5,df6)

  return(Final.df)
}

########################################################################################
# Expected length and sum of length using 6 base methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)
#
##### 10. Expected length for a given n and alpha level for 6 base methods
explAll<-function(n,alp,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")

  #### Calling functions and creating df
  df.1    = gexplWD(n,alp,a,b)
  df.2    = gexplSC(n,alp,a,b)
  df.3    = gexplAS(n,alp,a,b)
  df.4    = gexplLT(n,alp,a,b)
  df.5    = gexplTW(n,alp,a,b)
  df.6    = gexplLR(n,alp,a,b)

  df.new=  rbind(df.1,df.2,df.3,df.4,df.5,df.6)
  return(df.new)
}
