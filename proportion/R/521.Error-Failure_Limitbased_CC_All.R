#' Calculates error, long term power and pass/fail criteria for continuity corrected Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correction
#' @param phi - Null hypothesis value
#' @param f - Failure criterion
#' @details  Evaluation of Wald-type interval with continuity correction using error
#' due to the difference of achieved and nominal level of significance for the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{delalp}{ Delta-alpha is the increase of the nominal error with respect to real error}
#'  \item{theta}{ Long term power of the test}
#'  \item{Fail_Pass}{Fail/pass based on the input f criterion}
#' @family Error for continuity corrected methods
#' @examples
#' n=5; alp=0.05; phi=0.05;c=1/(2*n); f=-2
#' errCWD(n,alp,phi,c,f)
#' @references
#' [1] 2014 Martin Andres, A. and Alvarez Hernandez, M.
#' Two-tailed asymptotic inferences for a proportion.
#' Journal of Applied Statistics, 41, 7, 1516-1529
#' @export
##### 1.CC WALD - DELTA_ALPHA, THETA,F-ERROR,POWER,FAILURE
errCWD<-function(n,alp,phi,c,f)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(phi)) stop("'phi' is missing")
  if (missing(c)) stop("'c' is missing")
  if (missing(f)) stop("'f' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if (phi>1 || phi<0) stop("Null hypothesis 'phi' has to be between 0 and 1")
  if ((class(c) != "integer") & (class(c) != "numeric") || length(c) >1 || c<0 ) stop("'c' has to be positive")
  if ((class(f) != "integer") & (class(f) != "numeric")) stop("'f' has to be numeric value")

####INPUT n
x=0:n
k=n+1
####INITIALIZATIONS
pCW=0
qCW=0
seCW=0
LCW=0
UCW=0

###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)

for(i in 1:k)
{
pCW[i]=x[i]/n
qCW[i]=1-pCW[i]
seCW[i]=sqrt(pCW[i]*qCW[i]/n)
LCW[i]=max(pCW[i]-((cv*seCW[i])+c),0)
UCW[i]=min(pCW[i]+((cv*seCW[i])+c),1)
}

###DELTA_ALPHA, THETA,F
#z_alp=(qnorm(1-(alp/2),0,1))^2
alpstarCW=0
thetactr=0
for(m in 1:k)
{
if(phi > UCW[m] || phi<LCW[m])
{
thetactr=thetactr+1
alpstarCW[m]=dbinom(x[m],n,phi)
} else alpstarCW[m] = 0
}
delalpCW=round((alp-sum(alpstarCW))*100,2)
theta=round(100*thetactr/(n+1),2)
if(delalpCW<f)
Fail_Pass="failure" else Fail_Pass="success"
data.frame(delalp=delalpCW,theta,Fail_Pass)
}
#####################################################################################################################################
#' Calculates error, long term power and pass/fail criteria for continuity corrected Score method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correction
#' @param phi - Null hypothesis value
#' @param f - Failure criterion
#' @details  Evaluation of continuity corrected score test approach using error due to the
#' difference of achieved and nominal level of significance for the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{delalp}{ Delta-alpha is the increase of the nominal error with respect to real error}
#'  \item{theta}{ Long term power of the test}
#'  \item{Fail_Pass}{Fail/pass based on the input f criterion}
#' @family Error for continuity corrected methods
#' @examples
#' n=5; alp=0.05; phi=0.05;c=1/(2*n); f=-2
#' errCSC(n,alp,phi,c,f)
#' @references
#' [1] 2014 Martin Andres, A. and Alvarez Hernandez, M.
#' Two-tailed asymptotic inferences for a proportion.
#' Journal of Applied Statistics, 41, 7, 1516-1529
#' @export
##### 2.CC-SCORE:DELTA_ALPHA, THETA,F-ERROR,POWER,FAILURE
errCSC<-function(n,alp,phi,c,f)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(phi)) stop("'phi' is missing")
  if (missing(c)) stop("'c' is missing")
  if (missing(f)) stop("'f' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if (phi>1 || phi<0) stop("Null hypothesis 'phi' has to be between 0 and 1")
  if (c<=0 || c>(1/(2*n))) stop("'c' has to be positive and less than or equal to 1/(2*n)")
  if ((class(f) != "integer") & (class(f) != "numeric")) stop("'f' has to be numeric value")

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
LCS[i]=max((n/(n+(cv)^2))*((pCS[i]-c+cv1)-(cv2*seCS_L[i])),0)
UCS[i]=min((n/(n+(cv)^2))*((pCS[i]+c+cv1)+(cv2*seCS_U[i])),1)
}

###DELTA_ALPHA, THETA,F
#z_alp=(qnorm(1-(alp/2),0,1))^2
alpstarCS=0
thetactr=0
for(m1 in 1:k)
{
if(phi > UCS[m1] || phi<LCS[m1])
{
thetactr=thetactr+1
alpstarCS[m1]=dbinom(x[m1],n,phi)
} else alpstarCS[m1] = 0
}

delalpCS=round((alp-sum(alpstarCS))*100,2)
theta=round(100*thetactr/(n+1),2)
if(delalpCS<f)
Fail_Pass="failure" else Fail_Pass="success"
return(data.frame(delalp=delalpCS,theta,Fail_Pass))
}

#####################################################################################################################################
#' Calculates error, long term power and pass/fail criteria for continuity corrected ArcSine method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correction
#' @param phi - Null hypothesis value
#' @param f - Failure criterion
#' @details  Evaluation of continuity corrected Wald-type interval for the arcsine
#' transformation of the parameter \code{p}
#' using error due to the difference of achieved and nominal level of significance for the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{delalp}{ Delta-alpha is the increase of the nominal error with respect to real error}
#'  \item{theta}{ Long term power of the test}
#'  \item{Fail_Pass}{Fail/pass based on the input f criterion}
#' @family Error for continuity corrected methods
#' @examples
#' n=5; alp=0.05; phi=0.05;c=1/(2*n); f=-2
#' errCAS(n,alp,phi,c,f)
#' @references
#' [1] 2014 Martin Andres, A. and Alvarez Hernandez, M.
#' Two-tailed asymptotic inferences for a proportion.
#' Journal of Applied Statistics, 41, 7, 1516-1529
#' @export
##### 3.CC ARC SINE:DELTA_ALPHA, THETA,F-ERROR,POWER,FAILURE
errCAS<-function(n,alp,phi,c,f)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(phi)) stop("'phi' is missing")
  if (missing(c)) stop("'c' is missing")
  if (missing(f)) stop("'f' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if (phi>1 || phi<0) stop("Null hypothesis 'phi' has to be between 0 and 1")
  if ((class(c) != "integer") & (class(c) != "numeric") || length(c) >1 || c<0 ) stop("'c' has to be positive")
  if ((class(f) != "integer") & (class(f) != "numeric")) stop("'f' has to be numeric value")

####INPUT n
x=0:n
k=n+1
####INITIALIZATIONS
pCA=0
qCA=0
seCA=0
LCA=0
UCA=0
###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)

#ARC-SINE METHOD
for(i in 1:k)
{
pCA[i]=x[i]/n
qCA[i]=1-pCA[i]
seCA[i]=cv/sqrt(4*n)
LCA[i]=max((sin(asin(sqrt(pCA[i]))-seCA[i]-c))^2,0)
UCA[i]=min((sin(asin(sqrt(pCA[i]))+seCA[i]+c))^2,1)
}

###DELTA_ALPHA, THETA,F
#z_alp=(qnorm(1-(alp/2),0,1))^2
alpstarCA=0
thetactr=0
for(m1 in 1:k)
{
if(phi > UCA[m1] || phi<LCA[m1])
{
thetactr=thetactr+1
alpstarCA[m1]=dbinom(x[m1],n,phi)
} else alpstarCA[m1] = 0
}

delalpCA=round((alp-sum(alpstarCA))*100,2)
theta=round(100*thetactr/(n+1),2)
if(delalpCA<f)
Fail_Pass="failure" else Fail_Pass="success"
return(data.frame(delalp=delalpCA,theta,Fail_Pass))
}

#####################################################################################################################################
#' Calculates error, long term power and pass/fail criteria for continuity corrected Logit Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correction
#' @param phi - Null hypothesis value
#' @param f - Failure criterion
#' @details  Evaluation of continuity corrected Wald-type interval based on the logit transformation of \code{p}
#' using error due to the difference of achieved and nominal level of significance for the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{delalp}{ Delta-alpha is the increase of the nominal error with respect to real error}
#'  \item{theta}{ Long term power of the test}
#'  \item{Fail_Pass}{Fail/pass based on the input f criterion}
#' @family Error for continuity corrected methods
#' @examples
#' n=5; alp=0.05; phi=0.05;c=1/(2*n); f=-2
#' errCLT(n,alp,phi,c,f)
#' @references
#' [1] 2014 Martin Andres, A. and Alvarez Hernandez, M.
#' Two-tailed asymptotic inferences for a proportion.
#' Journal of Applied Statistics, 41, 7, 1516-1529
#' @export
##### 4.CC LOGIT:DEALTA_ALPHA, THETA,F-ERROR,POWER,FAILURE
errCLT<-function(n,alp,phi,c,f)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(phi)) stop("'phi' is missing")
  if (missing(c)) stop("'c' is missing")
  if (missing(f)) stop("'f' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if (phi>1 || phi<0) stop("Null hypothesis 'phi' has to be between 0 and 1")
  if ((class(c) != "integer") & (class(c) != "numeric") || length(c) >1 || c<0 ) stop("'c' has to be positive")
  if ((class(f) != "integer") & (class(f) != "numeric")) stop("'f' has to be numeric value")


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
LCLT[j+1]=max(lgiti(lgit[j+1]-(cv/seCLT[j+1])-c),0)
UCLT[j+1]=min(lgiti(lgit[j+1]+(cv/seCLT[j+1])+c),1)
}
#z_alp=(qnorm(1-(alp/2),0,1))^2
alpstarCLT=0
thetactr=0
for(m in 1:k)
{
if(phi > UCLT[m] || phi<LCLT[m])
{
thetactr=thetactr+1
alpstarCLT[m]=dbinom(x[m],n,phi)
} else alpstarCLT[m] = 0
}

delalpCLT=round((alp-sum(alpstarCLT))*100,2)
theta=round(100*thetactr/(n+1),2)
if(delalpCLT<f)
Fail_Pass="failure" else Fail_Pass="success"
return(data.frame(delalp=delalpCLT,theta,Fail_Pass))
}

#####################################################################################################################################
#' Calculates error, long term power and pass/fail criteria for continuity corrected Wald-T method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correction
#' @param phi - Null hypothesis value
#' @param f - Failure criterion
#' @details  Evaluation of approximate and continuity corrected method based on a t_approximation of the standardized point estimator
#' using error due to the difference of achieved and nominal level of significance for the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{delalp}{ Delta-alpha is the increase of the nominal error with respect to real error}
#'  \item{theta}{ Long term power of the test}
#'  \item{Fail_Pass}{Fail/pass based on the input f criterion}
#' @family Error for continuity corrected methods
#' @examples
#' n=5; alp=0.05; phi=0.05;c=1/(2*n); f=-2
#' errCTW(n,alp,phi,c,f)
#' @references
#' [1] 2014 Martin Andres, A. and Alvarez Hernandez, M.
#' Two-tailed asymptotic inferences for a proportion.
#' Journal of Applied Statistics, 41, 7, 1516-1529
#' @export
##### 5.CC WALD_t:DELTA_ALPHA, THETA,F-ERROR,POWER,FAILURE
errCTW<-function(n,alp,phi,c,f)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(phi)) stop("'phi' is missing")
  if (missing(c)) stop("'c' is missing")
  if (missing(f)) stop("'f' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if (phi>1 || phi<0) stop("Null hypothesis 'phi' has to be between 0 and 1")
  if ((class(c) != "integer") & (class(c) != "numeric") || length(c) >1 || c<0 ) stop("'c' has to be positive")
  if ((class(f) != "integer") & (class(f) != "numeric")) stop("'f' has to be numeric value")

####DATA
x=0:n
k=n+1
####INITIALIZATIONS
pCTW=0
qCTW=0
seCTW=0
DOF=0
cv=0
LCTW=0
UCTW=0
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
LCTW[i]=max(pCTW[i]-(seCTW[i]+c),0)
UCTW[i]=min(pCTW[i]+(seCTW[i]+c),1)
}
#z_alp=(qnorm(1-(alp/2),0,1))^2
alpstarCTW=0
thetactr=0
for(m in 1:k)
{
if(phi > UCTW[m] || phi<LCTW[m])
{
thetactr=thetactr+1
alpstarCTW[m]=dbinom(x[m],n,phi)
} else alpstarCTW[m] = 0
}

delalpCTW=round((alp-sum(alpstarCTW))*100,2)
theta=round(100*thetactr/(n+1),2)
if(delalpCTW<f)
Fail_Pass="failure" else Fail_Pass="success"
return(data.frame(delalp=delalpCTW,theta,Fail_Pass))
}

#########################################################################################
#' Calculates error, long term power and pass/fail criteria using 5 continuity corrected methods (Wald, Wald-T, Score, Logit-Wald, ArcSine)
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continuity correction
#' @param phi - Null hypothesis value
#' @param f - Failure criterion
#' @details  Calculates error, long term power and pass/fail criteria using 5
#' continuity corrected methods (Wald, Wald-T, Score, Logit-Wald, ArcSine)
#' @return A dataframe with
#'  \item{delalp}{ Delta-alpha is the increase of the nominal error with respect to real error}
#'  \item{theta}{ Long term power of the test}
#'  \item{Fail_Pass}{Fail/pass based on the input f criterion}
#'  \item{method}{Name of the method}
#' @family Error for continuity corrected methods
#' @examples
#' n=5; alp=0.05; phi=0.05;c=1/(2*n); f=-2
#' errCAll(n,alp,phi,c,f)
#' @references
#' [1] 2014 Martin Andres, A. and Alvarez Hernandez, M.
#' Two-tailed asymptotic inferences for a proportion.
#' Journal of Applied Statistics, 41, 7, 1516-1529
#' @export
##### 1.CC WALD - DELTA_ALPHA, THETA,F-ERROR,POWER,FAILURE
errCAll<-function(n,alp,phi,c,f)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(phi)) stop("'phi' is missing")
  if (missing(c)) stop("'c' is missing")
  if (missing(f)) stop("'f' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if (phi>1 || phi<0) stop("Null hypothesis 'phi' has to be between 0 and 1")
  if (c<=0 || c>(1/(2*n))) stop("'c' has to be positive and less than or equal to 1/(2*n)")
  if ((class(f) != "integer") & (class(f) != "numeric")) stop("'f' has to be numeric value")

  df.1    = errCWD(n,alp,phi,c,f)
  df.2    = errCSC(n,alp,phi,c,f)
  df.3    = errCAS(n,alp,phi,c,f)
  df.4    = errCLT(n,alp,phi,c,f)
  df.5    = errCTW(n,alp,phi,c,f)

  df.1$method = as.factor("CC-Wald")
  df.2$method = as.factor("CC-Score")
  df.3$method = as.factor("CC-ArcSine")
  df.4$method = as.factor("CC-Logit-Wald")
  df.5$method = as.factor("CC-Wald-T")

  df.new=  rbind(df.1,df.2,df.3,df.4,df.5)
  return(df.new)

}
