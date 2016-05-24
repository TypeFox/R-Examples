#' Calculates error, long term power and pass/fail criteria for Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param phi - Null hypothesis value
#' @param f - Failure criterion
#' @details  Evaluation of Wald-type intervals using error due to the
#' difference of achieved and nominal level of significance for the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{delalp}{ Delta-alpha is the increase of the nominal error with respect to real error}
#'  \item{theta}{ Long term power of the test}
#'  \item{Fail_Pass}{Fail/pass based on the input f criterion}
#' @family Error for base methods
#' @examples
#' n=20; alp=0.05; phi=0.05; f=-2
#' errWD(n,alp,phi,f)
#' @references
#' [1] 2014 Martin Andres, A. and Alvarez Hernandez, M.
#' Two-tailed asymptotic inferences for a proportion.
#' Journal of Applied Statistics, 41, 7, 1516-1529
#' @export
##### 1.WALD - DELTA_ALPHA, THETA,F-ERROR,POWER,FAILURE
errWD<-function(n,alp,phi,f)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(phi)) stop("'phi' is missing")
  if (missing(f)) stop("'f' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if (phi>1 || phi<0 || length(phi)>1) stop("Null hypothesis 'phi' has to be between 0 and 1")
  if ((class(f) != "integer") & (class(f) != "numeric")|| length(f)>1) stop("'f' has to be numeric value")

####DATA
x=0:n
k=n+1
####INITIALIZATIONS
pW=0
qW=0
seW=0
LWD=0
UWD=0
alpstarW=0
thetactr=0

###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#WALD METHOD
for(i in 1:k)
{
pW[i]=x[i]/n
qW[i]=1-(x[i]/n)
seW[i]=sqrt(pW[i]*qW[i]/n)
LWD[i]=max(pW[i]-(cv*seW[i]),0)
UWD[i]=min(pW[i]+(cv*seW[i]),1)
}
for(m in 1:k)
{
if(phi > UWD[m] || phi<LWD[m])
{
thetactr=thetactr+1
alpstarW[m]=dbinom(x[m],n,phi)
} else alpstarW[m] = 0
}
delalpW=round((alp-sum(alpstarW))*100,2)
theta=round(100*thetactr/(n+1),2)
if(delalpW<f)
Fail_Pass="failure" else Fail_Pass="success"
data.frame(delalp=delalpW,theta,Fail_Pass)
}

#####################################################################################################################################
#' Calculates error, long term power and pass/fail criteria for Score method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param phi - Null hypothesis value
#' @param f - Failure criterion
#' @details  Evaluation of score test approach using error due to the
#' difference of achieved and nominal level of significance for the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{delalp}{ Delta-alpha is the increase of the nominal error with respect to real error}
#'  \item{theta}{ Long term power of the test}
#'  \item{Fail_Pass}{Fail/pass based on the input f criterion}
#' @family Error for base methods
#' @examples
#' n=20; alp=0.05; phi=0.05; f=-2
#' errSC(n,alp,phi,f)
#' @references
#' [1] 2014 Martin Andres, A. and Alvarez Hernandez, M.
#' Two-tailed asymptotic inferences for a proportion.
#' Journal of Applied Statistics, 41, 7, 1516-1529
#' @export
##### 2.SCORE:DELTA_ALPHA, THETA,F-ERROR,POWER,FAILURE
errSC<-function(n,alp,phi,f)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(phi)) stop("'phi' is missing")
  if (missing(f)) stop("'f' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if (phi>1 || phi<0 || length(phi)>1) stop("Null hypothesis 'phi' has to be between 0 and 1")
  if ((class(f) != "integer") & (class(f) != "numeric")|| length(f)>1) stop("'f' has to be numeric value")

####DATA
x=0:n
k=n+1
####INITIALIZATIONS
pS=0
qS=0
seS=0
LSC=0
USC=0

#SCORE (WILSON) METHOD
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
LSC[i]=max((n/(n+(cv)^2))*((pS[i]+cv1)-(cv*seS[i])),0)
USC[i]=min((n/(n+(cv)^2))*((pS[i]+cv1)+(cv*seS[i])),1)
}
###DELTA_ALPHA, THETA,F
alpstarS=0
thetactr=0
for(m in 1:k)
{
if(phi > USC[m] || phi<LSC[m])
{
thetactr=thetactr+1
alpstarS[m]=dbinom(x[m],n,phi)
} else alpstarS[m] = 0
}

delalpS=round((alp-sum(alpstarS))*100,2)
theta=round(100*thetactr/(n+1),2)
if(delalpS<f)
Fail_Pass="failure" else Fail_Pass="success"
return(data.frame(delalp=delalpS,theta,Fail_Pass))
}
#####################################################################################################################################
#' Calculates error, long term power and pass/fail criteria for ArcSine method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param phi - Null hypothesis value
#' @param f - Failure criterion
#' @details  Evaluation of Wald-type interval for the arcsine transformation of the parameter
#' \code{p} error due to the difference of achieved and nominal level of
#' significance for the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{delalp}{ Delta-alpha is the increase of the nominal error with respect to real error}
#'  \item{theta}{ Long term power of the test}
#'  \item{Fail_Pass}{Fail/pass based on the input f criterion}
#' @family Error for base methods
#' @examples
#' n=20; alp=0.05; phi=0.05; f=-2
#' errAS(n,alp,phi,f)
#' @references
#' [1] 2014 Martin Andres, A. and Alvarez Hernandez, M.
#' Two-tailed asymptotic inferences for a proportion.
#' Journal of Applied Statistics, 41, 7, 1516-1529
#' @export
##### 3.ARC SINE:DELTA_ALPHA, THETA,F-ERROR,POWER,FAILURE
errAS<-function(n,alp,phi,f)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(phi)) stop("'phi' is missing")
  if (missing(f)) stop("'f' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if (phi>1 || phi<0 || length(phi)>1) stop("Null hypothesis 'phi' has to be between 0 and 1")
  if ((class(f) != "integer") & (class(f) != "numeric")|| length(f)>1) stop("'f' has to be numeric value")

####DATA
x=0:n
k=n+1
####INITIALIZATIONS
pA=0
qA=0
seA=0
LAS=0
UAS=0
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#ARC-SINE METHOD
for(i in 1:k)
{
pA[i]=x[i]/n
qA[i]=1-pA[i]
seA[i]=cv/sqrt(4*n)
LAS[i]=max((sin(asin(sqrt(pA[i]))-seA[i]))^2,0)
UAS[i]=min((sin(asin(sqrt(pA[i]))+seA[i]))^2,1)
}
###DELTA_ALPHA, THETA,F
alpstarA=0
thetactr=0
for(m in 1:k)
{
if(phi > UAS[m] || phi<LAS[m])
{
thetactr=thetactr+1
alpstarA[m]=dbinom(x[m],n,phi)
} else alpstarA[m] = 0
}
delalpA=round((alp-sum(alpstarA))*100,2)
theta=round(100*thetactr/(n+1),2)
if(delalpA<f)
Fail_Pass="failure" else Fail_Pass="success"
return(data.frame(delalp=delalpA,theta,Fail_Pass))
}

#####################################################################################################################################
#' Calculates error, long term power and pass/fail criteria for Logit Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param phi - Null hypothesis value
#' @param f - Failure criterion
#' @details  Evaluation of Wald-type interval based on the logit transformation of \code{p}
#' using error due to the difference of achieved and nominal level of significance for the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{delalp}{ Delta-alpha is the increase of the nominal error with respect to real error}
#'  \item{theta}{ Long term power of the test}
#'  \item{Fail_Pass}{Fail/pass based on the input f criterion}
#' @family Error for base methods
#' @examples
#' n=20; alp=0.05; phi=0.05; f=-2
#' errLT(n,alp,phi,f)
#' @references
#' [1] 2014 Martin Andres, A. and Alvarez Hernandez, M.
#' Two-tailed asymptotic inferences for a proportion.
#' Journal of Applied Statistics, 41, 7, 1516-1529
#' @export
##### 4.LOGIT WALD :DELTA_ALPHA, THETA,F-ERROR,POWER,FAILURE
errLT<-function(n,alp,phi,f)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(phi)) stop("'phi' is missing")
  if (missing(f)) stop("'f' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if (phi>1 || phi<0 || length(phi)>1) stop("Null hypothesis 'phi' has to be between 0 and 1")
  if ((class(f) != "integer") & (class(f) != "numeric")|| length(f)>1) stop("'f' has to be numeric value")

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
}
###DELTA_ALPHA, THETA,F
alpstarLT=0
thetactr=0
for(m in 1:k)
{
if(phi > ULT[m] || phi<LLT[m])
{
thetactr=thetactr+1
alpstarLT[m]=dbinom(x[m],n,phi)
} else alpstarLT[m] = 0
}

delalpLT=round((alp-sum(alpstarLT))*100,2)
theta=round(100*thetactr/(n+1),2)
if(delalpLT<f)
Fail_Pass="failure" else Fail_Pass="success"
return(data.frame(delalp=delalpLT,theta,Fail_Pass))
}

#####################################################################################################################################
#' Calculates error, long term power and pass/fail criteria for Wald-T method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param phi - Null hypothesis value
#' @param f - Failure criterion
#' @details  Evaluation of approximate method based on a t_approximation of the
#' standardized point estimator using error due to the difference of achieved and nominal level of significance for the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{delalp}{ Delta-alpha is the increase of the nominal error with respect to real error}
#'  \item{theta}{ Long term power of the test}
#'  \item{Fail_Pass}{Fail/pass based on the input f criterion}
#' @family Error for base methods
#' @examples
#' n=20; alp=0.05; phi=0.05; f=-2
#' errTW(n,alp,phi,f)
#' @references
#' [1] 2014 Martin Andres, A. and Alvarez Hernandez, M.
#' Two-tailed asymptotic inferences for a proportion.
#' Journal of Applied Statistics, 41, 7, 1516-1529
#' @export
##### 5. WALD -t:DELTA_ALPHA, THETA,F-ERROR,POWER,FAILURE
errTW<-function(n,alp,phi,f)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(phi)) stop("'phi' is missing")
  if (missing(f)) stop("'f' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if (phi>1 || phi<0 || length(phi)>1) stop("Null hypothesis 'phi' has to be between 0 and 1")
  if ((class(f) != "integer") & (class(f) != "numeric")|| length(f)>1) stop("'f' has to be numeric value")

####DATA
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
LTW[i]=max(pTW[i]-(seTW[i]),0)
UTW[i]=min(pTW[i]+(seTW[i]),1)
}
###DELTA_ALPHA, THETA,F
alpstarTW=0
thetactr=0
for(m in 1:k)
{
if(phi > UTW[m] || phi<LTW[m])
{
thetactr=thetactr+1
alpstarTW[m]=dbinom(x[m],n,phi)
} else alpstarTW[m] = 0
}

delalpTW=round((alp-sum(alpstarTW))*100,2)
theta=round(100*thetactr/(n+1),2)
if(delalpTW<f)
Fail_Pass="failure" else Fail_Pass="success"
return(data.frame(delalp=delalpTW,theta,Fail_Pass))
}

#####################################################################################################################################
#' Calculates error, long term power and pass/fail criteria for Likelihood Ratio method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param phi - Null hypothesis value
#' @param f - Failure criterion
#' @details  Evaluation of Likelihood ratio limits using error due to the difference
#' of achieved and nominal level of significance for the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{delalp}{ Delta-alpha is the increase of the nominal error with respect to real error}
#'  \item{theta}{ Long term power of the test}
#'  \item{Fail_Pass}{Fail/pass based on the input f criterion}
#' @family Error for base methods
#' @examples
#' n=20; alp=0.05; phi=0.05; f=-2
#' errLR(n,alp,phi,f)
#' @references
#' [1] 2014 Martin Andres, A. and Alvarez Hernandez, M.
#' Two-tailed asymptotic inferences for a proportion.
#' Journal of Applied Statistics, 41, 7, 1516-1529
#' @export
##### 6.LIKELIHOOD RATIO:DELTA_ALPHA, THETA,F-ERROR,POWER,FAILURE
errLR<-function(n,alp,phi,f)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(phi)) stop("'phi' is missing")
  if (missing(f)) stop("'f' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if (phi>1 || phi<0 || length(phi)>1) stop("Null hypothesis 'phi' has to be between 0 and 1")
  if ((class(f) != "integer") & (class(f) != "numeric")|| length(f)>1) stop("'f' has to be numeric value")

####DATA
y=0:n
k=n+1
####INITIALIZATIONS
mle=0
cutoff=0
LLR=0
ULR=0

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
LLR[i]=optimize(loglik.optim, c(0,mle[i]))$minimum
ULR[i]=optimize(loglik.optim, c(mle[i],1))$minimum
}
###DELTA_ALPHA, THETA,F
alpstarL=0
thetactr=0
for(m in 1:k)
{
if(phi > ULR[m] || phi < LLR[m])
{
thetactr=thetactr+1
alpstarL[m]=dbinom(y[m],n,phi)
} else alpstarL[m] = 0
}

delalpL=round((alp-sum(alpstarL))*100,2)
theta=round(100*thetactr/(n+1),2)
if(delalpL<f)
Fail_Pass="failure" else Fail_Pass="success"
return(data.frame(delalp=delalpL,theta,Fail_Pass))
}

#####################################################################################################################################
#' Calculates error, long term power and pass/fail criteria for Exact method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param phi - Null hypothesis value
#' @param f - Failure criterion
#' @param e - Exact method indicator  in [0, 1] {1: Clopper Pearson, 0.5: Mid P}
#' The input can also be a range of values between 0 and 1.
#' @details  Evaluation of Confidence interval for \code{p}
#' based on inverting equal-tailed binomial tests with null hypothesis \eqn{H0: p = p0}
#' using error due to the difference of achieved and nominal level of significance for the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{delalp}{ Delta-alpha is the increase of the nominal error with respect to real error}
#'  \item{theta}{ Long term power of the test}
#'  \item{Fail_Pass}{Fail/pass based on the input f criterion}
#' @family Error for base methods
#' @examples
#' n=20; alp=0.05;phi=0.05; f=-2;e=0.5 # Mid-p
#' errEX(n,alp,phi,f,e)
#' n=20; alp=0.05;phi=0.05; f=-2;e=1 #Clopper-Pearson
#' errEX(n,alp,phi,f,e)
#' n=20; alp=0.05;phi=0.05; f=-2;e=c(0.1,0.5,0.95,1) #Range including Mid-p and Clopper-Pearson
#' errEX(n,alp,phi,f,e)
#' @references
#' [1] 2014 Martin Andres, A. and Alvarez Hernandez, M.
#' Two-tailed asymptotic inferences for a proportion.
#' Journal of Applied Statistics, 41, 7, 1516-1529
#' @export
#7.EXACT METHOD
errEX<-function(n,alp,phi,f,e)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(phi)) stop("'phi' is missing")
  if (missing(f)) stop("'f' is missing")
  if (missing(e)) stop("'e' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if (phi>1 || phi<0 || length(phi)>1) stop("Null hypothesis 'phi' has to be between 0 and 1")
  if ((class(f) != "integer") & (class(f) != "numeric")|| length(f)>1) stop("'f' has to be numeric value")
  if ((class(e) != "integer") & (class(e) != "numeric") || any(e>1) || any(e<0)) stop("'e' has to be between 0 and 1")
  if (length(e)>10) stop("'e' can have only 10 intervals")

  nvar=length(e)

  res <- data.frame()

  for(i in 1:nvar)
  {
    lu=gerrEX501(n,alp,phi,f,e[i])
    res <- rbind(res,lu)
  }
  return(res)
}
gerrEX501<-function(n,alp,phi,f,e)
{
x=0:n
k=n+1
#####Exact Limits
LEX=0
UEX=0
for(i in 1:k)
{
LEX[i]=exlim501l(x[i],n,alp,e)
UEX[i]=exlim501u(x[i],n,alp,e)
}

###DELTA_ALPHA, THETA,F
alpstarE=0
thetactr=0
for(m in 1:k)
{
if(phi > UEX[m] || phi < LEX[m])
{
thetactr=thetactr+1
alpstarE[m]=dbinom(x[m],n,phi)
} else alpstarE[m] = 0
}

delalpE=round((alp-sum(alpstarE))*100,2)
theta=round(100*thetactr/(n+1),2)
if(delalpE<f)
Fail_Pass="failure" else Fail_Pass="success"
return(data.frame(delalp=delalpE,theta,Fail_Pass,e))
}

exlim501l=function(x,n,alp,e)
{
  if(x==0)
  {
    LEX = 0
  } else if(x==n){
    LEX= (alp/(2*e))^(1/n)
  }else
  {
    z=x-1
    y=0:z
    f1=function(p) (1-e)*dbinom(x,n,p)+sum(dbinom(y,n,p))-(1-(alp/2))
    LEX= uniroot(f1,c(0,1))$root
  }
  return(LEX)
}
exlim501u=function(x,n,alp,e)
{
  if(x==0)
  {
    UEX = 1-((alp/(2*e))^(1/n))
  } else if(x==n){
    UEX = 1
  }else
  {
    z=x-1
    y=0:z
    f2= function(p) e*dbinom(x,n,p)+sum(dbinom(y,n,p))-(alp/2)
    UEX =uniroot(f2,c(0,1))$root
  }
  return(UEX)
}


######################################################################
#' Calculates error, long term power and pass/fail criteria for Bayesian method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param phi - Null hypothesis value
#' @param f - Failure criterion
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  Evaluation of Bayesian Highest Probability Density
#' (HPD) and two tailed intervals using error due to the difference of achieved and
#' nominal level of significance for the \eqn{n + 1} intervals
#' for the Beta - Binomial conjugate prior model for the probability of success \code{p}
#' @return A dataframe with
#'  \item{delalp}{ Delta-alpha is the increase of the nominal error with respect to real error}
#'  \item{theta}{ Long term power of the test}
#'  \item{Fail_Pass}{Fail/pass based on the input f criterion}
#'  \item{method}{Name of method - Quantile or HPD}
#' @family Error for base methods
#' @examples
#' n=20; alp=0.05; phi=0.05; f=-2;a=0.5;b=0.5
#' errBA(n,alp,phi,f,a,b)
#' @references
#' [1] 2014 Martin Andres, A. and Alvarez Hernandez, M.
#' Two-tailed asymptotic inferences for a proportion.
#' Journal of Applied Statistics, 41, 7, 1516-1529
#' @export
#8.BAYESIAN
errBA<-function(n,alp,phi,f,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(phi)) stop("'phi' is missing")
  if (missing(f)) stop("'f' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if (phi>1 || phi<0 || length(phi)>1) stop("Null hypothesis 'phi' has to be between 0 and 1")
  if ((class(f) != "integer") & (class(f) != "numeric")|| length(f)>1) stop("'f' has to be numeric value")
  if ((class(a) != "integer") & (class(a) != "numeric") || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || b<0  ) stop("'b' has to be greater than or equal to 0")

####INPUT n
x=0:n
k=n+1
####INITIALIZATIONS
LBAQ=0
UBAQ=0
LBAH=0
UBAH=0
##############
#library(TeachingDemos)				#To get HPDs
for(i in 1:k)
{
#Quantile Based Intervals
LBAQ[i]=qbeta(alp/2,x[i]+a,n-x[i]+b)
UBAQ[i]=qbeta(1-(alp/2),x[i]+a,n-x[i]+b)

LBAH[i]=TeachingDemos::hpd(qbeta,shape1=x[i]+a,shape2=n-x[i]+b,conf=1-alp)[1]
UBAH[i]=TeachingDemos::hpd(qbeta,shape1=x[i]+a,shape2=n-x[i]+b,conf=1-alp)[2]

}
###DELTA_ALPHA, THETA,F_Quantile Based
alpstarLQ=0
thetactrQ=0
for(m in 1:k)
{
if(phi > UBAQ[m] || phi < LBAQ[m])
{
thetactrQ=thetactrQ+1
alpstarLQ[m]=dbinom(x[m],n,phi)
} else alpstarLQ[m] = 0
}

delalpLQ=round((alp-sum(alpstarLQ))*100,2)
thetaQ=round(100*thetactrQ/(n+1),2)
if(delalpLQ<f)
Fail_PassQ="failure" else Fail_PassQ="success"

###DELTA_ALPHA, THETA,F_HPD Based
alpstarLH=0
thetactrH=0
for(m in 1:k)
{
if(phi > UBAH[m] || phi < LBAH[m])
{
thetactrH=thetactrH+1
alpstarLH[m]=dbinom(x[m],n,phi)
} else alpstarLH[m] = 0
}

delalpLH=round((alp-sum(alpstarLH))*100,2)
thetaH=round(100*thetactrH/(n+1),2)
if(delalpLH<f)
Fail_PassH="failure" else Fail_PassH="success"
qdf=data.frame(delalp=delalpLQ,theta=thetaQ,Fail_Pass=Fail_PassQ,method="Quantile")
hdf=data.frame(delalp=delalpLH,theta=thetaH,Fail_Pass=Fail_PassH,method="HPD")
ndf=rbind(qdf,hdf)
return(ndf)
}

########################################################################################
#' Calculates error, long term power and pass/fail criteria using 6 base methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param phi - Null hypothesis value
#' @param f - Failure criterion
#' @details  Calculation of error, long term power and pass/fail
#' criteria using 6 base methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)
#' @return A dataframe with
#'  \item{delalp}{ Delta-alpha is the increase of the nominal error with respect to real error}
#'  \item{theta}{ Long term power of the test}
#'  \item{Fail_Pass}{Fail/pass based on the input f criterion}
#'  \item{method}{Name of the method }
#' @family Error for base methods
#' @examples
#' n=20; alp=0.05; phi=0.05; f=-2
#' errAll(n,alp,phi,f)
#' @references
#' [1] 2014 Martin Andres, A. and Alvarez Hernandez, M.
#' Two-tailed asymptotic inferences for a proportion.
#' Journal of Applied Statistics, 41, 7, 1516-1529
#' @export
##### 10. Expected length for a given n and alpha level for 6 base methods
errAll<-function(n,alp,phi,f)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(phi)) stop("'phi' is missing")
  if (missing(f)) stop("'f' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if (phi>1 || phi<0 || length(phi)>1) stop("Null hypothesis 'phi' has to be between 0 and 1")
  if ((class(f) != "integer") & (class(f) != "numeric")|| length(f)>1) stop("'f' has to be numeric value")

  #### Calling functions and creating df
  df.1    = errWD(n,alp,phi,f)
  df.2    = errSC(n,alp,phi,f)
  df.3    = errAS(n,alp,phi,f)
  df.4    = errLT(n,alp,phi,f)
  df.5    = errTW(n,alp,phi,f)
  df.6    = errLR(n,alp,phi,f)

  df.1$method = as.factor("Wald")
  df.2$method = as.factor("Score")
  df.3$method = as.factor("ArcSine")
  df.4$method = as.factor("Logit-Wald")
  df.5$method = as.factor("Wald-T")
  df.6$method = as.factor("Likelihood")

  df.new=  rbind(df.1,df.2,df.3,df.4,df.5,df.6)
  return(df.new)
}
