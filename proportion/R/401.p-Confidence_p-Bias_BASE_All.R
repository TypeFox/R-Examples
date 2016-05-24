#' p-confidence and p-bias for Wald method given n and alpha level
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @details  Evaluation of Wald-type intervals using p-confidence and p-bias for the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{x1}{  Number of successes (positive samples)}
#'  \item{pconf }{   p-Confidence}
#'  \item{pbias }{   p-Bias}
#' @family p-confidence and p-bias of base methods
#' @examples
#' n=5; alp=0.05
#' pCOpBIWD(n,alp)
#' @references
#' [1] 2005 Vos PW and Hudson S.
#' Evaluation Criteria for Discrete Confidence Intervals: Beyond Coverage and Length.
#' The American Statistician: 59; 137 - 142.
#' @export
##### 1.WALD- p-confidence and p-bias for a given n and alpha level
pCOpBIWD<-function(n,alp)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")

  ####INPUT n
  x=0:n
  k=n+1
  ####INITIALIZATIONS
  pW=0
  qW=0
  seW=0
  LW=0
  UW=0
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
    pW[i]=x[i]/n
    qW[i]=1-(x[i]/n)
    seW[i]=sqrt(pW[i]*qW[i]/n)
    LW[i]=pW[i]-(cv*seW[i])
    UW[i]=pW[i]+(cv*seW[i])
    if(LW[i]<0) LW[i]=0
    if(UW[i]>1) UW[i]=1
  }
  ####p-confidence and p-bias
  for(i in 2:(k-1))
  {
    pcon[i-1]=2*(pbinom(i-1, n, LW[i], lower.tail = FALSE, log.p = FALSE)+dbinom(i-1, n, LW[i]))
    pconC[i-1]=2*pbinom(i-1, n, UW[i], lower.tail = TRUE, log.p = FALSE)
    pconf[i-1]=(1-max(pcon[i-1],pconC[i-1]))*100 		#p-confidence calculation
    pbia1[i-1]=max(pcon[i-1],pconC[i-1])-min(pcon[i-1],pconC[i-1])
    pbias[i-1]=max(0,pbia1[i-1])*100
  }
  x1=1:(n-1)
  p_C_B=data.frame(x1,pconf,pbias)
#   windows()
#   plot(x1,pconf,type="l",xlab="No of successes",ylab="p-Confidence",main="Wald")
#   windows()
#   plot(x1,pbias,type="l",xlab="No of successes",ylab="p-Bias",main="Wald")
  return(p_C_B)
}
#####################################################################################################################################
#' p-confidence and p-bias for  Score method given n and alpha level
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @details  Evaluation of score test approach using p-confidence and p-bias for the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{x1}{  Number of successes (positive samples)}
#'  \item{pconf }{   p-Confidence}
#'  \item{pbias }{   p-Bias}
#' @family p-confidence and p-bias of base methods
#' @examples
#' n=5; alp=0.05
#' pCOpBISC(n,alp)
#' @references
#' [1] 2005 Vos PW and Hudson S.
#' Evaluation Criteria for Discrete Confidence Intervals: Beyond Coverage and Length.
#' The American Statistician: 59; 137 - 142.
#' @export
##### 2.SCORE- p-confidence and p-bias for a given n and alpha level
pCOpBISC<-function(n,alp)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")

  ####INPUT n
  x=0:n
  k=n+1
  ####INITIALIZATIONS
  pS=0
  qS=0
  seS=0
  LS=0
  US=0
  pcon=0						#p-confidence
  pconC=0
  pconf=0
  pbia1=0					#p-bias
  pbias=0
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
  }
  ####p-confidence and p-bias
  for(i in 2:(k-1))
  {
    pcon[i-1]=2*(pbinom(i-1, n, LS[i], lower.tail = FALSE, log.p = FALSE)+dbinom(i-1, n, LS[i]))
    pconC[i-1]=2*pbinom(i-1, n, US[i], lower.tail = TRUE, log.p = FALSE)
    pconf[i-1]=1-max(pcon[i-1],pconC[i-1]) 		#p-confidence calculation
    pbia1[i-1]=max(pcon[i-1],pconC[i-1])-min(pcon[i-1],pconC[i-1])
    pbias[i-1]=as.numeric(max(0,pbia1[i-1]))
  }
  x1=1:(n-1)
  p_C_B=data.frame(x1,pconf,pbias)
#   windows()
#   plot(x1,pconf,type="l",xlab="No of successes",ylab="p-Confidence",main="Score")
#   windows()
#   plot(x1,pbias,type="l",xlab="No of successes",ylab="p-Bias",main="Score")
  return(p_C_B)
}
#####################################################################################################################################
#' p-confidence and p-bias for  ArcSine method given n and alpha level
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @details  Evaluation of Wald-type interval for the arcsine transformation of the parameter \code{p} using p-confidence and p-bias for the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{x1}{  Number of successes (positive samples)}
#'  \item{pconf }{   p-Confidence}
#'  \item{pbias }{   p-Bias}
#' @family p-confidence and p-bias of base methods
#' @examples
#' n=5; alp=0.05
#' pCOpBIAS(n,alp)
#' @references
#' [1] 2005 Vos PW and Hudson S.
#' Evaluation Criteria for Discrete Confidence Intervals: Beyond Coverage and Length.
#' The American Statistician: 59; 137 - 142.
#' @export
##### 3.ARC SINE - p-confidence and p-bias for a given n and alpha level
pCOpBIAS<-function(n,alp)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")

  ####INPUT n
  x=0:n
  k=n+1
  ####INITIALIZATIONS
  pA=0
  qA=0
  seA=0
  LA=0
  UA=0
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
    pA[i]=x[i]/n
    qA[i]=1-pA[i]
    seA[i]=cv/sqrt(4*n)
    LA[i]=(sin(asin(sqrt(pA[i]))-seA[i]))^2
    UA[i]=(sin(asin(sqrt(pA[i]))+seA[i]))^2
    if(LA[i]<0) LA[i]=0
    if(UA[i]>1) UA[i]=1
  }

  ####p-confidence and p-bias
  for(i in 2:(k-1))
  {
    pcon[i-1]=2*(pbinom(i-1, n, LA[i], lower.tail = FALSE, log.p = FALSE)+dbinom(i-1, n, LA[i]))
    pconC[i-1]=2*pbinom(i-1, n, UA[i], lower.tail = TRUE, log.p = FALSE)
    pconf[i-1]=1-max(pcon[i-1],pconC[i-1]) 		#p-confidence calculation
    pbia1[i-1]=max(pcon[i-1],pconC[i-1])-min(pcon[i-1],pconC[i-1])
    pbias[i-1]=as.numeric(max(0,pbia1[i-1]))
  }

  x1=1:(n-1)
  p_C_B=data.frame(x1,pconf,pbias)
  # windows()
  # plot(x1,pconf,type="l",xlab="No of successes",ylab="p-Confidence",main="Arc Sine")
  # windows()
  # plot(x1,pbias,type="l",xlab="No of successes",ylab="p-Bias",main="Arc Sine")
  return(p_C_B)
}
#####################################################################################################################################
#' p-confidence and p-bias for  Logit Wald method given n and alpha level
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @details  Evaluation of Wald-type interval based on the logit transformation of \code{p} using p-confidence and p-bias for the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{x1}{  Number of successes (positive samples)}
#'  \item{pconf }{   p-Confidence}
#'  \item{pbias }{   p-Bias}
#' @family p-confidence and p-bias of base methods
#' @examples
#' n=5; alp=0.05
#' pCOpBILT(n,alp)
#' @references
#' [1] 2005 Vos PW and Hudson S.
#' Evaluation Criteria for Discrete Confidence Intervals: Beyond Coverage and Length.
#' The American Statistician: 59; 137 - 142.
#' @export
##### 4.LOGIT WALD - p-confidence and p-bias for a given n and alpha level
pCOpBILT<-function(n,alp)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")

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
  pcon=0						#p-confidence
  pconC=0
  pconf=0
  pbia1=0					#p-bias
  pbias=0
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


  ####p-confidence and p-bias
  for(i in 2:(k-1))
  {
    pcon[i-1]=2*(pbinom(i-1, n, LLT[i], lower.tail = FALSE, log.p = FALSE)+dbinom(i-1, n, LLT[i]))
    pconC[i-1]=2*pbinom(i-1, n, ULT[i], lower.tail = TRUE, log.p = FALSE)
    pconf[i-1]=1-max(pcon[i-1],pconC[i-1]) 		#p-confidence calculation
    pbia1[i-1]=max(pcon[i-1],pconC[i-1])-min(pcon[i-1],pconC[i-1])
    pbias[i-1]=as.numeric(max(0,pbia1[i-1]))
  }
  x1=1:(n-1)
  p_C_B=data.frame(x1,pconf,pbias)
  # windows()
  # plot(x1,pconf,type="l",xlab="No of successes",ylab="p-Confidence",main="Logit Wald")
  # windows()
  # plot(x1,pbias,type="l",xlab="No of successes",ylab="p-Bias",main="Logit Wald")
  return(p_C_B)
}
#####################################################################################################################################
#' p-confidence and p-bias for  T-Wald method given n and alpha level
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @details  Evaluation of approximate method based on a t_approximation of the standardized point estimator using p-confidence and p-bias for the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{x1}{  Number of successes (positive samples)}
#'  \item{pconf }{   p-Confidence}
#'  \item{pbias }{   p-Bias}
#' @family p-confidence and p-bias of base methods
#' @examples
#' n=5; alp=0.05
#' pCOpBITW(n,alp)
#' @references
#' [1] 2005 Vos PW and Hudson S.
#' Evaluation Criteria for Discrete Confidence Intervals: Beyond Coverage and Length.
#' The American Statistician: 59; 137 - 142.
#' @export
##### 5.T-WALD: p-confidence and p-bias for a given n and alpha level
pCOpBITW<-function(n,alp)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")

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
  }
  ####p-confidence and p-bias
  for(i in 2:(k-1))
  {
    pcon[i-1]=2*(pbinom(i-1, n, LTW[i], lower.tail = FALSE, log.p = FALSE)+dbinom(i-1, n, LTW[i]))
    pconC[i-1]=2*pbinom(i-1, n, UTW[i], lower.tail = TRUE, log.p = FALSE)
    pconf[i-1]=1-max(pcon[i-1],pconC[i-1]) 		#p-confidence calculation
    pbia1[i-1]=max(pcon[i-1],pconC[i-1])-min(pcon[i-1],pconC[i-1])
    pbias[i-1]=as.numeric(max(0,pbia1[i-1]))
  }
  x1=1:(n-1)
  p_C_B=data.frame(x1,pconf,pbias)
  # windows()
  # plot(x1,pconf,type="l",xlab="No of successes",ylab="p-Confidence",main="t-Wald")
  # windows()
  # plot(x1,pbias,type="l",,xlab="No of successes",ylab="p-Bias",main="t-Wald")
  return(p_C_B)
}

#####################################################################################################################################
#' p-confidence and p-bias for Likelihood method given n and alpha level
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @details  Evaluation of Likelihood ratio limits using p-confidence and p-bias for the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{x1}{  Number of successes (positive samples)}
#'  \item{pconf }{   p-Confidence}
#'  \item{pbias }{   p-Bias}
#' @family p-confidence and p-bias of base methods
#' @examples
#' n=5; alp=0.05
#' pCOpBILR(n,alp)
#' @references
#' [1] 2005 Vos PW and Hudson S.
#' Evaluation Criteria for Discrete Confidence Intervals: Beyond Coverage and Length.
#' The American Statistician: 59; 137 - 142.
#' @export
##### 6.LIKELIHOOD RATIO - p-confidence and p-bias for a given n and alpha level
pCOpBILR<-function(n,alp)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")

  ####INPUT n
  y=0:n
  k=n+1
  ####INITIALIZATIONS
  mle=0
  cutoff=0
  LL=0
  UL=0
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
    likelhd = function(p) dbinom(y[i],n,p)
    loglik = function(p) dbinom(y[i],n,p,log=TRUE)
    mle[i]=optimize(likelhd,c(0,1),maximum=TRUE)$maximum
    cutoff[i]=loglik(mle[i])-(cv^2/2)
    loglik.optim=function(p){abs(cutoff[i]-loglik(p))}
    LL[i]=optimize(loglik.optim, c(0,mle[i]))$minimum
    UL[i]=optimize(loglik.optim, c(mle[i],1))$minimum
  }
  ####p-confidence and p-bias
  for(i in 2:(k-1))
  {
    pcon[i-1]=2*(pbinom(i-1, n, LL[i], lower.tail = FALSE, log.p = FALSE)+dbinom(i-1, n, LL[i]))
    pconC[i-1]=2*pbinom(i-1, n, UL[i], lower.tail = TRUE, log.p = FALSE)
    pconf[i-1]=1-max(pcon[i-1],pconC[i-1]) 		#p-confidence calculation
    pbia1[i-1]=max(pcon[i-1],pconC[i-1])-min(pcon[i-1],pconC[i-1])
    pbias[i-1]=as.numeric(max(0,pbia1[i-1]))
  }
  x1=1:(n-1)
  p_C_B=data.frame(x1,pconf,pbias)
  # windows()
  # plot(x1,pconf,type="l",xlab="No of successes",ylab="p-Confidence",main="Likelihood Ratio")
  # windows()
  # plot(x1,pbias,type="l",xlab="No of successes",ylab="p-Bias",main="Likelihood Ratio")
  return(p_C_B)
}
#####################################################################################################################################
#' p-confidence and p-bias for  Exact method given n and alpha level
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param e - Exact method indicator  in [0, 1] {1: Clopper Pearson, 0.5: Mid P}.
#' The input can also be a range of values between 0 and 1.
#' @details  Evaluation of Confidence interval for \code{p} based on inverting equal-tailed binomial tests with null hypothesis \eqn{H0: p = p0} using p-confidence and p-bias for the \eqn{n + 1} intervals
#' @return A dataframe with
#'  \item{x1}{  Number of successes (positive samples)}
#'  \item{pconf }{   p-Confidence}
#'  \item{pbias }{   p-Bias}
#'  \item{e }{- Exact method input}
#' @family p-confidence and p-bias of base methods
#' @examples
#' n=5; alp=0.05;e=0.5
#' pCOpBIEX(n,alp,e)
#' n=5; alp=0.05;e=1 #Clopper-Pearson
#' pCOpBIEX(n,alp,e)
#' n=5; alp=0.05;e=c(0.1,0.5,0.95,1) #Range including Mid-p and Clopper-Pearson
#' pCOpBIEX(n,alp,e)
#' @references
#' [1] 2005 Vos PW and Hudson S.
#' Evaluation Criteria for Discrete Confidence Intervals: Beyond Coverage and Length.
#' The American Statistician: 59; 137 - 142.
#' @export
##### 7.EXACT METHODS - p-confidence and p-bias for a given n and alpha level
pCOpBIEX<-function(n,alp,e)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(e)) stop("'e' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(e) != "integer") & (class(e) != "numeric") || any(e>1) || any(e<0)) stop("'e' has to be between 0 and 1")
  if (length(e)>10) stop("'e' can have only 10 intervals")

  ####INITIALIZATIONS
  nvar=length(e)

  tdf <- data.frame()

  for(i in 1:nvar)
     {
     slf=slu401(n,alp,e[i])
     tdf <- rbind(tdf,slf)
     }

  res <- data.frame()
  stdf <- data.frame()

  for(i in 1:nvar)
      {
       ep=e[i]
       stdf=subset(tdf, e == ep)
       pcbw=pcb401(n,stdf)
       res <- rbind(res,pcbw)
      }
    return(res)
}

slu401<-function(n,alp,e)
{
  x=0:n
  k=n+1
  LEX=0
  UEX=0
  #EXACT METHOD
  LEX[1]=0
  UEX[k]=1
  for(i in 1:(k-2))
  {
    UEX[1]= 1-((alp/(2*e))^(1/n))
    LEX[k]=(alp/(2*e))^(1/n)
    LEX[i+1]=exlim401l(x[i+1],n,alp,e)
    UEX[i+1]=exlim401u(x[i+1],n,alp,e)

    }
  return(data.frame(LEX,UEX,e))
}

pcb401=function(n,stdf)
{
#  x=0:n
  k=n+1
  pcon=0						#p-confidence
  pconC=0
  pconf=0
  pbia1=0					#p-bias
  pbias=0
  LEX=stdf$LEX
  UEX=stdf$UEX
  e=stdf$e

  for(i in 2:(k-1))
  {
    pcon[i-1]=2*(pbinom(i-1, n, LEX[i], lower.tail = FALSE, log.p = FALSE)+dbinom(i-1, n, LEX[i]))
    pconC[i-1]=2*pbinom(i-1, n, UEX[i], lower.tail = TRUE, log.p = FALSE)
    pconf[i-1]=1-max(pcon[i-1],pconC[i-1]) 		#p-confidence calculation
    pbia1[i-1]=max(pcon[i-1],pconC[i-1])-min(pcon[i-1],pconC[i-1])
    pbias[i-1]=as.numeric(max(0,pbia1[i-1]))
  }
  x1=1:(n-1)
  p_C_B=data.frame(x1,pconf,pbias,e=e[1])
  return(p_C_B)

}
#####TO FIND LOWER LIMITS
exlim401l=function(x,n,alp,e)
{
  z=x-1
  y=0:z
  f1=function(p) (1-e)*dbinom(x,n,p)+sum(dbinom(y,n,p))-(1-(alp/2))
  LEX= uniroot(f1,c(0,1))$root
  return(LEX)
}
#####TO FIND UPPER LIMITS
exlim401u=function(x,n,alp,e)
{
  z=x-1
  y=0:z
  f2  = function(p) e*dbinom(x,n,p)+sum(dbinom(y,n,p))-(alp/2)
  UEX = uniroot(f2,c(0,1))$root
  return(UEX)
}

#######################################################################################################
#' p-confidence and p-bias for  Bayesian method given n and alpha level and priors a & b
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a1 - Shape parameter 1 for prior Beta distribution in Bayesian model
#' @param a2 - Shape parameter 2 for prior Beta distribution in Bayesian model
#' @details  Evaluation of Bayesian Highest Probability Density (HPD) and two tailed
#' intervals using p-confidence and p-bias for the \eqn{n + 1} intervals for the
#' Beta - Binomial conjugate prior model for the probability of success \code{p}
#' @return A dataframe with
#'  \item{x1}{  Number of successes (positive samples)}
#'  \item{pconfQ }{   p-Confidence Quantile}
#'  \item{pbiasQ }{   p-Bias Quantile}
#'  \item{pconfH }{   p-Confidence HPD}
#'  \item{pbiasH }{   p-Bias HPD}
#' @family p-confidence and p-bias of base methods
#' @examples
#' n=5; alp=0.05;a1=1;a2=1
#' pCOpBIBA(n,alp,a1,a2)
#' @references
#' [1] 2005 Vos PW and Hudson S.
#' Evaluation Criteria for Discrete Confidence Intervals: Beyond Coverage and Length.
#' The American Statistician: 59; 137 - 142.
#' @export
####8.BAYESIAN p-confidence and p-bias
pCOpBIBA<-function(n,alp,a1,a2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a1)) stop("'a1' is missing")
  if (missing(a2)) stop("'a2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
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
  pconQ=0
  pconCQ=0
  pconfQ=0
  pbia1Q=0
  pbiasQ=0
  pconH=0
  pconCH=0
  pconfH=0
  pbia1H=0
  pbiasH=0

  ##############
  #library(TeachingDemos)				#To get HPDs
  for(i in 1:k)
  {
    #Quantile Based Intervals
    LBAQ[i]=qbeta(alp/2,x[i]+a1,n-x[i]+a2)
    UBAQ[i]=qbeta(1-(alp/2),x[i]+a1,n-x[i]+a2)

    LBAH[i]=TeachingDemos::hpd(qbeta,shape1=x[i]+a1,shape2=n-x[i]+a2,conf=1-alp)[1]
    UBAH[i]=TeachingDemos::hpd(qbeta,shape1=x[i]+a1,shape2=n-x[i]+a2,conf=1-alp)[2]
  }

  for(i in 2:(k-1))
  {
    pconQ[i-1]=2*(pbinom(i-1, n, LBAQ[i], lower.tail = FALSE, log.p = FALSE)+dbinom(i-1, n, LBAQ[i]))
    pconCQ[i-1]=2*pbinom(i-1, n, UBAQ[i], lower.tail = TRUE, log.p = FALSE)
    pconfQ[i-1]=(1-max(pconQ[i-1],pconCQ[i-1]))*100 		#p-confidence calculation
    pbia1Q[i-1]=max(pconQ[i-1],pconCQ[i-1])-min(pconQ[i-1],pconCQ[i-1])
    pbiasQ[i-1]=max(0,pbia1Q[i-1])*100

    pconH[i-1]=2*(pbinom(i-1, n, LBAH[i], lower.tail = FALSE, log.p = FALSE)+dbinom(i-1, n, LBAH[i]))
    pconCH[i-1]=2*pbinom(i-1, n, UBAH[i], lower.tail = TRUE, log.p = FALSE)
    pconfH[i-1]=(1-max(pconH[i-1],pconCH[i-1]))*100 		#p-confidence calculation
    pbia1H[i-1]=max(pconH[i-1],pconCH[i-1])-min(pconH[i-1],pconCH[i-1])
    pbiasH[i-1]=max(0,pbia1H[i-1])*100

  }
  x1=1:(n-1)

  return(data.frame(x1,pconfQ,pbiasQ,pconfH,pbiasH))
}

#####################################################################################
#' Calculates p-confidence and p-bias for a given n and alpha level for 6 base methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @details  Evaluation of  p-confidence and p-bias for  the \eqn{n + 1} intervals using 6 base methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)
#' @return A dataframe with
#'  \item{x1}{  Number of successes (positive samples)}
#'  \item{pconf }{   p-Confidence}
#'  \item{pbias }{   p-Bias}
#'  \item{method}{   Method}
#' @family p-confidence and p-bias of base methods
#' @examples
#' n=5; alp=0.05
#' pCOpBIAll(n,alp)
#' @references
#' [1] 2005 Vos PW and Hudson S.
#' Evaluation Criteria for Discrete Confidence Intervals: Beyond Coverage and Length.
#' The American Statistician: 59; 137 - 142.
#' @export
#10.All methods
pCOpBIAll<-function(n,alp)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")

  #### Calling functions and creating df

  df1 = pCOpBIWD(n,alp)
  df2 = pCOpBISC(n,alp)
  df3 = pCOpBIAS(n,alp)
  df4 = pCOpBILT(n,alp)
  df5 = pCOpBITW(n,alp)
  df6 = pCOpBILR(n,alp)

  df1$method = as.factor("Wald")
  df2$method = as.factor("Score")
  df3$method = as.factor("ArcSine")
  df4$method = as.factor("Logit-Wald")
  df5$method = as.factor("Wald-T")
  df6$method = as.factor("Likelihood")

  Final.df= rbind(df1,df2,df3,df4,df5,df6)

  return(Final.df)
}

