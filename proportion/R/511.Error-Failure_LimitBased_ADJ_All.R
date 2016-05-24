#' Calculates error, long term power and pass/fail criteria for adjusted Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param phi - Null hypothesis value
#' @param f - Failure criterion
#' @details  Evaluation of adjusted Wald-type interval using error due to the
#' difference of achieved and nominal level of significance for the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{delalp}{ Delta-alpha is the increase of the nominal error with respect to real error}
#'  \item{theta}{ Long term power of the test}
#'  \item{Fail_Pass}{Fail/pass based on the input f criterion}
#' @family Error for adjusted methods
#' @examples
#' n=20; alp=0.05; h=2;phi=0.99; f=-2
#' errAWD(n,alp,h,phi,f)
#' @references
#' [1] 2014 Martin Andres, A. and Alvarez Hernandez, M.
#' Two-tailed asymptotic inferences for a proportion.
#' Journal of Applied Statistics, 41, 7, 1516-1529
#' @export
##### 1.ADJUSTED WALD - DELTA_ALPHA, THETA,F-ERROR,POWER,FAILURE
errAWD<-function(n,alp,h,phi,f)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(phi)) stop("'phi' is missing")
  if (missing(f)) stop("'f' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || h<0  ) stop("'h' has to be greater than or equal to 0")
  if (phi>1 || phi<0) stop("Null hypothesis 'phi' has to be between 0 and 1")
  if ((class(f) != "integer") & (class(f) != "numeric")) stop("'f' has to be numeric value")

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


###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#WALD METHOD
for(i in 1:k)
{
pAW[i]=y[i]/n1
qAW[i]=1-pAW[i]
seAW[i]=sqrt(pAW[i]*qAW[i]/n1)
LAWD[i]=max(pAW[i]-(cv*seAW[i]),0)
UAWD[i]=min(pAW[i]+(cv*seAW[i]),1)
}

#####Finding Error, Failure
alpstarAW=0
thetactr=0
for(m in 1:k)
{
if(phi > UAWD[m] || phi<LAWD[m])
{
thetactr=thetactr+1
alpstarAW[m]=dbinom(x[m],n,phi)
} else alpstarAW[m] = 0
}
delalpAW=round((alp-sum(alpstarAW))*100,2)
theta=round(100*thetactr/(n+1),2)
if(delalpAW < f)Fail_Pass="failure" else Fail_Pass="success"
return(data.frame(delalp=delalpAW,theta,Fail_Pass))
}

##############################################################################
#' Calculates error, long term power and pass/fail criteria for adjusted Score method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param phi - Null hypothesis value
#' @param f - Failure criterion
#' @details  Evaluation of adjusted score test approach using error due to the
#' difference of achieved and nominal level of significance for the  \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{delalp}{ Delta-alpha is the increase of the nominal error with respect to real error}
#'  \item{theta}{ Long term power of the test}
#'  \item{Fail_Pass}{Fail/pass based on the input f criterion}
#' @family Error for adjusted methods
#' @examples
#' n=20; alp=0.05; h=2;phi=0.99; f=-2
#' errASC(n,alp,h,phi,f)
#' @references
#' [1] 2014 Martin Andres, A. and Alvarez Hernandez, M.
#' Two-tailed asymptotic inferences for a proportion.
#' Journal of Applied Statistics, 41, 7, 1516-1529
#' @export
#2.ADJUSTED SCORE
errASC<-function(n,alp,h,phi,f)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(phi)) stop("'phi' is missing")
  if (missing(f)) stop("'f' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || h<0  ) stop("'h' has to be greater than or equal to 0")
  if (phi>1 || phi<0) stop("Null hypothesis 'phi' has to be between 0 and 1")
  if ((class(f) != "integer") & (class(f) != "numeric")) stop("'f' has to be numeric value")

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
LASC[i]=max((n1/(n1+(cv)^2))*((pAS[i]+cv1)-(cv*seAS[i])),0)
UASC[i]=min((n1/(n1+(cv)^2))*((pAS[i]+cv1)+(cv*seAS[i])),1)
}
#####Finding Error, Failure
alpstarAS=0
thetactr=0
for(m in 1:k)
{
if(phi > UASC[m] || phi<LASC[m])
{
thetactr=thetactr+1
alpstarAS[m]=dbinom(x[m],n,phi)
} else alpstarAS[m] = 0
}
delalpAS=round((alp-sum(alpstarAS))*100,2)
theta=round(100*thetactr/(n+1),2)
if(delalpAS < f)Fail_Pass="failure" else Fail_Pass="success"
return(data.frame(delalp=delalpAS,theta,Fail_Pass))
}

##############################################################################
#' Calculates error, long term power and pass/fail criteria for adjusted ArcSine method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param phi - Null hypothesis value
#' @param f - Failure criterion
#' @details  Evaluation of adjusted Wald-type interval for the arcsine transformation of the parameter \code{p}
#' using error due to the difference of achieved and nominal level of significance for the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{delalp}{ Delta-alpha is the increase of the nominal error with respect to real error}
#'  \item{theta}{ Long term power of the test}
#'  \item{Fail_Pass}{Fail/pass based on the input f criterion}
#' @family Error for adjusted methods
#' @examples
#' n=20; alp=0.05; h=2;phi=0.99; f=-2
#' errAAS(n,alp,h,phi,f)
#' @references
#' [1] 2014 Martin Andres, A. and Alvarez Hernandez, M.
#' Two-tailed asymptotic inferences for a proportion.
#' Journal of Applied Statistics, 41, 7, 1516-1529
#' @export
#3.ADJUSTED ARCSINE
errAAS<-function(n,alp,h,phi,f)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(phi)) stop("'phi' is missing")
  if (missing(f)) stop("'f' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || h<0  ) stop("'h' has to be greater than or equal to 0")
  if (phi>1 || phi<0) stop("Null hypothesis 'phi' has to be between 0 and 1")
  if ((class(f) != "integer") & (class(f) != "numeric")) stop("'f' has to be numeric value")

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

cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#ARC-SINE METHOD
for(i in 1:k)
{
pA[i]=y[i]/m
qA[i]=1-pA[i]
seA[i]=cv/sqrt(4*m)
LAAS[i]= max((sin(asin(sqrt(pA[i]))-seA[i]))^2,0)
UAAS[i]= min((sin(asin(sqrt(pA[i]))+seA[i]))^2,1)
}
#####Finding Error, Failure
alpstarAAS=0
thetactr=0
for(m in 1:k)
{
if(phi > UAAS[m] || phi<LAAS[m])
{
thetactr=thetactr+1
alpstarAAS[m]=dbinom(x[m],n,phi)
} else alpstarAAS[m] = 0
}
delalpAAS=round((alp-sum(alpstarAAS))*100,2)
theta=round(100*thetactr/(n+1),2)
if(delalpAAS < f)Fail_Pass="failure" else Fail_Pass="success"
return(data.frame(delalp=delalpAAS,theta,Fail_Pass))
}

###########################################################################################################
#' Calculates error, long term power and pass/fail criteria for adjusted Likelihood Ratio method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param phi - Null hypothesis value
#' @param f - Failure criterion
#' @details  Evaluation of adjusted Likelihood ratio limits using error due to the
#' difference of achieved and nominal level of significance for the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{delalp}{ Delta-alpha is the increase of the nominal error with respect to real error}
#'  \item{theta}{ Long term power of the test}
#'  \item{Fail_Pass}{Fail/pass based on the input f criterion}
#' @family Error for adjusted methods
#' @examples
#' n=20; alp=0.05; h=2;phi=0.99; f=-2
#' errALR(n,alp,h,phi,f)
#' @references
#' [1] 2014 Martin Andres, A. and Alvarez Hernandez, M.
#' Two-tailed asymptotic inferences for a proportion.
#' Journal of Applied Statistics, 41, 7, 1516-1529
#' @export
#4.ADJUSTED LIKELIHOOD RATIO
errALR<-function(n,alp,h,phi,f)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(phi)) stop("'phi' is missing")
  if (missing(f)) stop("'f' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h) >1|| h<0  || !(h%%1 ==0)) stop("'h' has to be an integer greater than or equal to 0")
  if (phi>1 || phi<0) stop("Null hypothesis 'phi' has to be between 0 and 1")
  if ((class(f) != "integer") & (class(f) != "numeric")) stop("'f' has to be numeric value")

####INPUT n
y=0:n
y1=y+h
k=n+1
n1=n+(2*h)
####INITIALIZATIONS
mle=0
cutoff=0
LALR=0
UALR=0

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
LALR[i]=max(optimize(loglik.optim, c(0,mle[i]))$minimum,0)
UALR[i]=min(optimize(loglik.optim, c(mle[i],1))$minimum,1)
}
#####Finding Error, Failure
alpstarALR=0
thetactr=0
for(m in 1:k)
{
if(phi > UALR[m] || phi<LALR[m])
{
thetactr=thetactr+1
alpstarALR[m]=dbinom(y[m],n,phi)
} else alpstarALR[m] = 0
}
delalpALR=round((alp-sum(alpstarALR))*100,2)
theta=round(100*thetactr/(n+1),2)
if(delalpALR < f)Fail_Pass="failure" else Fail_Pass="success"
return(data.frame(delalp=delalpALR,theta,Fail_Pass))
}

###########################################################################################################
#' Calculates error, long term power and pass/fail criteria for adjusted Wald-T method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param phi - Null hypothesis value
#' @param f - Failure criterion
#' @details  Evaluation of approximate and adjusted method based on a
#' t_approximation of the standardized point estimator
#' using error due to the difference of achieved and nominal level of significance for the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{delalp}{ Delta-alpha is the increase of the nominal error with respect to real error}
#'  \item{theta}{ Long term power of the test}
#'  \item{Fail_Pass}{Fail/pass based on the input f criterion}
#' @family Error for adjusted methods
#' @examples
#' n=20; alp=0.05; h=2;phi=0.99; f=-2
#' errATW(n,alp,h,phi,f)
#' @references
#' [1] 2014 Martin Andres, A. and Alvarez Hernandez, M.
#' Two-tailed asymptotic inferences for a proportion.
#' Journal of Applied Statistics, 41, 7, 1516-1529
#' @export
#5.ADJUSTED WALD-T
errATW<-function(n,alp,h,phi,f)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(phi)) stop("'phi' is missing")
  if (missing(f)) stop("'f' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || h<0  ) stop("'h' has to be greater than or equal to 0")
  if (phi>1 || phi<0) stop("Null hypothesis 'phi' has to be between 0 and 1")
  if ((class(f) != "integer") & (class(f) != "numeric")) stop("'f' has to be numeric value")

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
LATW[i]=max(pATW[i]-(seATW[i]),0)
UATW[i]=min(pATW[i]+(seATW[i]),1)
}
#####Finding Error, Failure
alpstarATW=0
thetactr=0
for(m in 1:k)
{
if(phi > UATW[m] || phi<LATW[m])
{
thetactr=thetactr+1
alpstarATW[m]=dbinom(x[m],n,phi)
} else alpstarATW[m] = 0
}
delalpATW=round((alp-sum(alpstarATW))*100,2)
theta=round(100*thetactr/(n+1),2)
if(delalpATW < f)Fail_Pass="failure" else Fail_Pass="success"
return(data.frame(delalp=delalpATW,theta,Fail_Pass))
}

###########################################################################################################
#' Calculates error, long term power and pass/fail criteria for adjusted Logit Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param phi - Null hypothesis value
#' @param f - Failure criterion
#' @details  Evaluation of adjusted Wald-type interval based on the logit transformation of \code{p}
#' using error due to the difference of achieved and nominal level of significance for the  \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{delalp}{ Delta-alpha is the increase of the nominal error with respect to real error}
#'  \item{theta}{ Long term power of the test}
#'  \item{Fail_Pass}{Fail/pass based on the input f criterion}
#' @family Error for adjusted methods
#' @examples
#' n=20; alp=0.05; h=2;phi=0.99; f=-2
#' errALT(n,alp,h,phi,f)
#' @references
#' [1] 2014 Martin Andres, A. and Alvarez Hernandez, M.
#' Two-tailed asymptotic inferences for a proportion.
#' Journal of Applied Statistics, 41, 7, 1516-1529
#' @export
#6.ADJUSTED LOGIT-WALD
errALT<-function(n,alp,h,phi,f)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(phi)) stop("'phi' is missing")
  if (missing(f)) stop("'f' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || h<0  ) stop("'h' has to be greater than or equal to 0")
  if (phi>1 || phi<0) stop("Null hypothesis 'phi' has to be between 0 and 1")
  if ((class(f) != "integer") & (class(f) != "numeric")) stop("'f' has to be numeric value")

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
###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#LOGIT-WALD METHOD
for(i in 1:k)
{
pALT[i]=y[i]/n1
qALT[i]=1-pALT[i]
lgit[i]=log(pALT[i]/qALT[i])
seALT[i]=sqrt(pALT[i]*qALT[i]*n1)
LALT[i]=max(1/(1+exp(-lgit[i]+(cv/seALT[i]))),1)
UALT[i]=min(1/(1+exp(-lgit[i]-(cv/seALT[i]))),1)
}
#####Finding Error, Failure
alpstarALT=0
thetactr=0
for(m in 1:k)
{
if(phi > UALT[m] || phi<LALT[m])
{
thetactr=thetactr+1
alpstarALT[m]=dbinom(x[m],n,phi)
} else alpstarALT[m] = 0
}
delalpALT=round((alp-sum(alpstarALT))*100,2)
theta=round(100*thetactr/(n+1),2)
if(delalpALT < f)Fail_Pass="failure" else Fail_Pass="success"
return(data.frame(delalp=delalpALT,theta,Fail_Pass))
}

###########################################################################################################
#' Calculates error, long term power and pass/fail criteria using 6 adjusted methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param phi - Null hypothesis value
#' @param f - Failure criterion
#' @details  Calculates  error, long term power and pass/fail
#' criteria using 6 adjusted methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)
#' @return A dataframe with
#'  \item{delalp}{ Delta-alpha is the increase of the nominal error with respect to real error}
#'  \item{theta}{ Long term power of the test}
#'  \item{Fail_Pass}{Fail/pass based on the input f criterion}
#'  \item{method}{Name of the method}
#' @family Error for adjusted methods
#' @examples
#' n=20; alp=0.05;h=2; phi=0.99; f=-2
#' errAAll(n,alp,h,phi,f)
#' @export
##### 10. Expected length for a given n and alpha level for 6 base methods
errAAll<-function(n,alp,h,phi,f)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(phi)) stop("'phi' is missing")
  if (missing(f)) stop("'f' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h) >1|| h<0  || !(h%%1 ==0)) stop("'h' has to be an integer greater than or equal to 0")
  if (phi>1 || phi<0) stop("Null hypothesis 'phi' has to be between 0 and 1")
  if ((class(f) != "integer") & (class(f) != "numeric")) stop("'f' has to be numeric value")
  #### Calling functions and creating df
  df.1    = errAWD(n,alp,h,phi,f)
  df.2    = errASC(n,alp,h,phi,f)
  df.3    = errAAS(n,alp,h,phi,f)
  df.4    = errALT(n,alp,h,phi,f)
  df.5    = errATW(n,alp,h,phi,f)
  df.6    = errALR(n,alp,h,phi,f)

  df.1$method = as.factor("Adj-Wald")
  df.2$method = as.factor("Adj-Score")
  df.3$method = as.factor("Adj-ArcSine")
  df.4$method = as.factor("Adj-Logit-Wald")
  df.5$method = as.factor("Adj-Wald-T")
  df.6$method = as.factor("Adj-Likelihood")

  df.new=  rbind(df.1,df.2,df.3,df.4,df.5,df.6)
  return(df.new)

}

