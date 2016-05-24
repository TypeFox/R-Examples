#' Coverage Probability of Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  Evaluation of Wald-type interval using coverage probability,
#' root mean square statistic, and the proportion of proportion lies within the
#' desired level of coverage
#' @return A dataframe with
#'  \item{mcpW}{  Wald Coverage Probability}
#'  \item{micpW }{   Wald minimum coverage probability}
#'  \item{RMSE_N }{   Root Mean Square Error from nominal size}
#'  \item{RMSE_M }{   Root Mean Square Error for Coverage Probability}
#'  \item{RMSE_MI }{  Root Mean Square Error for minimum coverage probability}
#'  \item{tol }{   Required tolerance for coverage probability}
#' @family Basic coverage probability methods
#' @examples
#' n= 10; alp=0.05; a=1;b=1; t1=0.93;t2=0.97
#' covpWD(n,alp,a,b,t1,t2)
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
#1.WALD
##### 1.WALD-Coverage Probability
covpWD<-function(n,alp,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")

####INPUT n
x=0:n
k=n+1
####INITIALIZATIONS
pW=0
qW=0
seW=0
LW=0
UW=0
s=5000								#Simulation run to generate hypothetical p
cpW=matrix(0,k,s)
ctW=matrix(0,k,s)							#Cover Pbty quantity in sum
cppW=0								#Coverage probabilty
RMSE_N1=0
RMSE_M1=0
RMSE_Mi1=0
ctr=0
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
####COVERAGE PROBABILITIES
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
if(hp[j] > LW[i] && hp[j] < UW[i])
{
cpW[i,j]=dbinom(i-1, n,hp[j])
ctW[i,j]=1
}
}
cppW[j]=sum(cpW[,j])
RMSE_N1[j]=(cppW[j]-(1-alp))^2			#Root mean Square from nominal size
if(t1<cppW[j]&&cppW[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined
}
#CPW=data.frame(hp,cppW)
mcpW=mean(cppW)
micpW=min(cppW)					#Mean Cov Prob
RMSE_N=sqrt(mean(RMSE_N1))

#Root mean Square from min and mean CoPr
for (j in 1:s)
{
RMSE_M1[j]=(cppW[j]-mcpW)^2
RMSE_Mi1[j]=(cppW[j]-micpW)^2
}
RMSE_M=sqrt(mean(RMSE_M1))
RMSE_MI=sqrt(mean(RMSE_Mi1))
tol=100*ctr/s
return(data.frame(mcpW,micpW,RMSE_N,RMSE_M,RMSE_MI,tol))
}

#############################################################################################################
#' Coverage Probability of Score method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  Evaluation of score test approach using coverage probability,
#' root mean square statistic, and the proportion of proportion lies within the
#' desired level of coverage
#' @return A dataframe with
#'  \item{mcpS}{  Score Coverage Probability}
#'  \item{micpS }{   Score minimum coverage probability}
#'  \item{RMSE_N }{   Root Mean Square Error from nominal size}
#'  \item{RMSE_M }{   Root Mean Square Error for Coverage Probability}
#'  \item{RMSE_MI }{  Root Mean Square Error for minimum coverage probability}
#'  \item{tol }{   Required tolerance for coverage probability}
#' @family Basic coverage probability methods
#' @examples
#' n= 10; alp=0.05; a=1;b=1; t1=0.93;t2=0.97
#' covpSC(n,alp,a,b,t1,t2)
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
##### 2.SCORE - Coverage Probability
covpSC<-function(n,alp,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")

  ###INPUT n
x=0:n
k=n+1
####INITIALIZATIONS
pS=0
qS=0
seS=0
LS=0
US=0
s=5000								#Simulation run to generate hypothetical p
cpS=matrix(0,k,s)
ctS=matrix(0,k,s)							#Cover Pbty quantity in sum
cppS=0								#Coverage probabilty
RMSE_N1=0
RMSE_M1=0
RMSE_Mi1=0
ctr=0

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
####COVERAGE PROBABILITIES
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
if(hp[j] > LS[i] && hp[j] < US[i])
{
cpS[i,j]=dbinom(i-1, n,hp[j])
ctS[i,j]=1
}
}
cppS[j]=sum(cpS[,j])						#Coverage Probability
RMSE_N1[j]=(cppS[j]-(1-alp))^2			#Root mean Square from nominal size
if(t1<cppS[j]&&cppS[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined
}
#CPS=data.frame(hp,cppS)
mcpS=mean(cppS)							#Mean Cov Prob
micpS=min(cppS)							#Mean Cov Prob
RMSE_N=sqrt(mean(RMSE_N1))

#Root mean Square from min and mean CoPr
for (j in 1:s)
{
RMSE_M1[j]=(cppS[j]-mcpS)^2
RMSE_Mi1[j]=(cppS[j]-micpS)^2
}
RMSE_M=sqrt(mean(RMSE_M1))
RMSE_MI=sqrt(mean(RMSE_Mi1))
tol=100*ctr/s
return(data.frame(mcpS,micpS,RMSE_N,RMSE_M,RMSE_MI,tol))
}


#############################################################################################################
#' Coverage Probability of ArcSine method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  Evaluation of Wald-type interval for the arcsine
#' transformation of the parameter p using coverage probability, root mean square
#' statistic, and the proportion of proportion lies within the desired level of coverage
#' @return A dataframe with
#'  \item{mcpA}{  ArcSine Coverage Probability}
#'  \item{micpA }{   ArcSine minimum coverage probability}
#'  \item{RMSE_N }{   Root Mean Square Error from nominal size}
#'  \item{RMSE_M }{   Root Mean Square Error for Coverage Probability}
#'  \item{RMSE_MI }{  Root Mean Square Error for minimum coverage probability}
#'  \item{tol }{   Required tolerance for coverage probability}
#' @family Basic coverage probability methods
#' @examples
#' n= 10; alp=0.05; a=1;b=1; t1=0.93;t2=0.97
#' covpAS(n,alp,a,b,t1,t2)
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
##### 3.ARC SINE - Coverage Probability
covpAS<-function(n,alp,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")

  ####INPUT n
x=0:n
k=n+1
####INITIALIZATIONS
pA=0
qA=0
seA=0
LA=0
UA=0
s=5000								#Simulation run to generate hypothetical p
cpA=matrix(0,k,s)
ctA=matrix(0,k,s)							#Cover Pbty quantity in sum
cppA=0								#Coverage probabilty
RMSE_N1=0
RMSE_M1=0
RMSE_Mi1=0
ctr=0

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
####COVERAGE PROBABILITIES
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
if(hp[j] > LA[i] && hp[j] < UA[i])
{
cpA[i,j]=dbinom(i-1, n,hp[j])
ctA[i,j]=1
}
}
cppA[j]=sum(cpA[,j])						#Coverage Probability
RMSE_N1[j]=(cppA[j]-(1-alp))^2			#Root mean Square from nominal size
if(t1<cppA[j]&&cppA[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined
}
#CPA=data.frame(hp,cppA)
mcpA=mean(cppA)							#Mean Cov Prob
micpA=min(cppA)
RMSE_N=sqrt(mean(RMSE_N1))

#Root mean Square from min and mean CoPr
for (j in 1:s)
{
RMSE_M1[j]=(cppA[j]-mcpA)^2
RMSE_Mi1[j]=(cppA[j]-micpA)^2
}
RMSE_M=sqrt(mean(RMSE_M1))
RMSE_MI=sqrt(mean(RMSE_Mi1))
tol=100*ctr/s
return(data.frame(mcpA,micpA,RMSE_N,RMSE_M,RMSE_MI,tol))
}
#############################################################################################################
#' Coverage Probability of Logit Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  Evaluation of Wald-type interval based on the logit transformation
#' of p using coverage probability, root mean square statistic, and the proportion
#' of proportion lies within the desired level of coverage
#' @return A dataframe with
#'  \item{mcpLT}{  Logit Wald Coverage Probability}
#'  \item{micpLT }{   Logit Wald minimum coverage probability}
#'  \item{RMSE_N }{   Root Mean Square Error from nominal size}
#'  \item{RMSE_M }{   Root Mean Square Error for Coverage Probability}
#'  \item{RMSE_MI }{  Root Mean Square Error for minimum coverage probability}
#'  \item{tol }{   Required tolerance for coverage probability}
#' @family Basic coverage probability methods
#' @examples
#' n= 10; alp=0.05; a=1;b=1; t1=0.93;t2=0.97
#' covpLT(n,alp,a,b,t1,t2)
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
##### 4.LOGIT-WALD - Coverage Probability
covpLT<-function(n,alp,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")

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
s=5000								#Simulation run to generate hypothetical p
cpLT=matrix(0,k,s)
ctLT=matrix(0,k,s)							#Cover Pbty quantity in sum
cppLT=0								#Coverage probabilty
RMSE_N1=0
RMSE_M1=0
RMSE_Mi1=0
ctr=0

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
####COVERAGE PROBABILITIES
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
if(hp[j] > LLT[i] && hp[j] < ULT[i])
{
cpLT[i,j]=dbinom(i-1, n,hp[j])
ctLT[i,j]=1
}
}
cppLT[j]=sum(cpLT[,j])				#Coverage Probability
RMSE_N1[j]=(cppLT[j]-(1-alp))^2			#Root mean Square from nominal size
if(t1<cppLT[j]&&cppLT[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined

}
#CPLT=data.frame(hp,cppLT)
mcpLT=mean(cppLT)
micpLT=min(cppLT)					#Mean Cov Prob
RMSE_N=sqrt(mean(RMSE_N1))

#Root mean Square from min and mean CoPr
for (j in 1:s)
{
RMSE_M1[j]=(cppLT[j]-mcpLT)^2
RMSE_Mi1[j]=(cppLT[j]-micpLT)^2
}
RMSE_M=sqrt(mean(RMSE_M1))
RMSE_MI=sqrt(mean(RMSE_Mi1))
tol=100*ctr/s
return(data.frame(mcpLT,micpLT,RMSE_N,RMSE_M,RMSE_MI,tol))
}
#############################################################################################################
#' Coverage Probability of Wald-T method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  Evaluation of approximate method based on a t_approximation of
#' the standardized point estimator using coverage probability, root mean square
#' statistic, and the proportion of proportion lies within the desired level of coverage
#' @return A dataframe with
#'  \item{mcpTW}{  Wald-T Coverage Probability}
#'  \item{micpTW }{   Wald-T minimum coverage probability}
#'  \item{RMSE_N }{   Root Mean Square Error from nominal size}
#'  \item{RMSE_M }{   Root Mean Square Error for Coverage Probability}
#'  \item{RMSE_MI }{  Root Mean Square Error for minimum coverage probability}
#'  \item{tol }{   Required tolerance for coverage probability}
#' @family Basic coverage probability methods
#' @examples
#' n= 10; alp=0.05; a=1;b=1; t1=0.93;t2=0.97
#' covpTW(n,alp,a,b,t1,t2)
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
##### 5. WALD_t - Coverage Probability
covpTW<-function(n,alp,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")

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
s=5000								#Simulation run to generate hypothetical p
cpTW=matrix(0,k,s)
ctTW=matrix(0,k,s)							#Cover Pbty quantity in sum
cppTW=0
RMSE_N1=0
RMSE_M1=0
RMSE_Mi1=0
ctr=0

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
####COVERAGE PROBABILITIES
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
if(hp[j] > LTW[i] && hp[j] < UTW[i])
{
cpTW[i,j]=dbinom(i-1, n,hp[j])
ctTW[i,j]=1
}
}
cppTW[j]=sum(cpTW[,j])						#Coverage Probability
RMSE_N1[j]=(cppTW[j]-(1-alp))^2			#Root mean Square from nominal size
if(t1<cppTW[j]&&cppTW[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined
}
#CPTW=data.frame(hp,cppTW)
mcpTW=mean(cppTW)
micpTW=min(cppTW)					#Mean Cov Prob
RMSE_N=sqrt(mean(RMSE_N1))

#Root mean Square from min and mean CoPr
for (j in 1:s)
{
RMSE_M1[j]=(cppTW[j]-mcpTW)^2
RMSE_Mi1[j]=(cppTW[j]-micpTW)^2
}
RMSE_M=sqrt(mean(RMSE_M1))
RMSE_MI=sqrt(mean(RMSE_Mi1))
tol=100*ctr/s
return(data.frame(mcpTW,micpTW,RMSE_N,RMSE_M,RMSE_MI,tol))
}

#############################################################################################################
#' Coverage Probability of  likelihood  method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  Evaluation of Likelihood ratio limits using coverage probability, root
#' mean square statistic, and the proportion of proportion lies within the desired
#' level of coverage
#' @return A dataframe with
#'  \item{mcpL}{   likelihood  Coverage Probability}
#'  \item{micpL }{    likelihood  minimum coverage probability}
#'  \item{RMSE_N }{   Root Mean Square Error from nominal size}
#'  \item{RMSE_M }{   Root Mean Square Error for Coverage Probability}
#'  \item{RMSE_MI }{  Root Mean Square Error for minimum coverage probability}
#'  \item{tol }{   Required tolerance for coverage probability}
#' @family Basic coverage probability methods
#' @examples
#' n= 10; alp=0.05; a=1;b=1; t1=0.93;t2=0.97
#' covpLR(n,alp,a,b,t1,t2)
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
##### 6.LIKELIHOOD RATIO - Coverage Probability
covpLR<-function(n,alp,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")

  ####INPUT n
y=0:n
k=n+1
####INITIALIZATIONS
mle=0
cutoff=0
LL=0
UL=0
s=5000								#Simulation run to generate hypothetical p
cpL=matrix(0,k,s)
ctL=matrix(0,k,s)							#Cover Pbty quantity in sum
cppL=0								#Coverage probabilty
RMSE_N1=0
RMSE_M1=0
RMSE_Mi1=0
ctr=0

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
####COVERAGE PROBABILITIES
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
if(hp[j] > LL[i] && hp[j] < UL[i])
{
cpL[i,j]=dbinom(i-1, n,hp[j])
ctL[i,j]=1
}
}
cppL[j]=sum(cpL[,j])						#Coverage Probability
RMSE_N1[j]=(cppL[j]-(1-alp))^2			#Root mean Square from nominal size
if(t1<cppL[j]&&cppL[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined
}
#CPL=data.frame(hp,cppL)
mcpL=mean(cppL)
micpL=min(cppL)					#Mean Cov Prob
RMSE_N=sqrt(mean(RMSE_N1))

#Root mean Square from min and mean CoPr
for (j in 1:s)
{
RMSE_M1[j]=(cppL[j]-mcpL)^2
RMSE_Mi1[j]=(cppL[j]-micpL)^2
}
RMSE_M=sqrt(mean(RMSE_M1))
RMSE_MI=sqrt(mean(RMSE_Mi1))
tol=100*ctr/s
return(data.frame(mcpL,micpL,RMSE_N,RMSE_M,RMSE_MI,tol))
}

#############################################################################################################
#' Coverage Probability of Exact  method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param e -  Exact method indicator (1:Clop-Pear,0.5:MID-p).
#' The input can also be a range of values between 0 and 1.
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  Evaluation of Confidence interval for p based on inverting equal-tailed
#' binomial tests with null hypothesis \eqn{H0: p = p0} using coverage probability,
#' root mean square statistic, and the proportion of proportion lies within the
#' desired level of coverage.
#' @return A dataframe with
#'  \item{mcpEX}{  Exact  Coverage Probability}
#'  \item{micpEX }{   Exact  minimum coverage probability}
#'  \item{RMSE_N }{   Root Mean Square Error from nominal size}
#'  \item{RMSE_M }{   Root Mean Square Error for Coverage Probability}
#'  \item{RMSE_MI }{  Root Mean Square Error for minimum coverage probability}
#'  \item{tol }{   Required tolerance for coverage probability}
#'  \item{e }{- Exact method input}
#' @family Basic coverage probability methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; e=0.5; a=1;b=1; t1=0.93;t2=0.97 # Mid-p
#' covpEX(n,alp,e,a,b,t1,t2)
#' n= 10; alp=0.05; e=1; a=1;b=1; t1=0.93;t2=0.97 #Clop-Pear
#' covpEX(n,alp,e,a,b,t1,t2)
#' n=5; alp=0.05;
#' e=c(0.1,0.5,0.95,1) #Range including Mid-p and Clopper-Pearson
#' a=1;b=1; t1=0.93;t2=0.97
#' covpEX(n,alp,e,a,b,t1,t2)
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
##### 7. EXACT METHOD - Coverage Probability
#EXACT METHOD
covpEX=function(n,alp,e,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(e)) stop("'e' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(e) != "integer") & (class(e) != "numeric") || any(e>1) || any(e<0)) stop("'e' has to be between 0 and 1")
  if (length(e)>10) stop("'e' can have only 10 intervals")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")

  nvar=length(e)

  res <- data.frame()

  for(i in 1:nvar)
  {
    lu=oldEX201(n,alp,e[i],a,b,t1,t2)
    res <- rbind(res,lu)
  }
  return(res)
}
oldEX201=function(n,alp,e,a,b,t1,t2)
{

x=0:n
k=n+1
LEX=0
UEX=0
s=5000					#Simulation run to generate hypothetical p
cpEX=matrix(0,k,s)
ctEX=matrix(0,k,s)			#Cover Pbty quantity in sum
cppEX=0
RMSE_N1=0
RMSE_M1=0
RMSE_Mi1=0
ctr=0
#EXACT METHOD
LEX[1]=0
UEX[1]= 1-((alp/(2*e))^(1/n))
LEX[k]=(alp/(2*e))^(1/n)
UEX[k]=1

for(i in 1:(k-2))
{
LEX[i+1]=exlim201l(x[i+1],n,alp,e)
UEX[i+1]=exlim201u(x[i+1],n,alp,e)
}
####COVERAGE PROBABILITIES
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
if(hp[j] > LEX[i] && hp[j] < UEX[i])
{
cpEX[i,j]=dbinom(i-1, n,hp[j])
ctEX[i,j]=1
}
}
cppEX[j]=sum(cpEX[,j])						#Coverage Probability
RMSE_N1[j]=(cppEX[j]-(1-alp))^2			#Root mean Square from nominal size
if(t1<cppEX[j]&&cppEX[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined
}
#CPEX=data.frame(hp,cppEX)
mcpEX=mean(cppEX)							#Mean Cov Prob
micpEX=min(cppEX)							#Mean Cov Prob
RMSE_N=sqrt(mean(RMSE_N1))

#Root mean Square from min and mean CoPr
for (j in 1:s)
{
RMSE_M1[j]=(cppEX[j]-mcpEX)^2
RMSE_Mi1[j]=(cppEX[j]-micpEX)^2
}
RMSE_M=sqrt(mean(RMSE_M1))
RMSE_MI=sqrt(mean(RMSE_Mi1))
tol=100*ctr/s
return(data.frame(mcpEX,micpEX,RMSE_N,RMSE_M,RMSE_MI,tol,e))
}
#####TO FIND LOWER LIMITS
exlim201l=function(x,n,alp,e)
{
  z=x-1
  y=0:z
  f1=function(p) (1-e)*dbinom(x,n,p)+sum(dbinom(y,n,p))-(1-(alp/2))
  LEX= uniroot(f1,c(0,1))$root
  return(LEX)
}
#####TO FIND UPPER LIMITS
exlim201u=function(x,n,alp,e)
{
  z=x-1
  y=0:z
  f2  = function(p) e*dbinom(x,n,p)+sum(dbinom(y,n,p))-(alp/2)
  UEX = uniroot(f2,c(0,1))$root
  return(UEX)
}


#############################################################################################################
#' Coverage Probability of Bayesian  method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @param a1 - Beta Prior Parameters for Bayesian estimation
#' @param a2 - Beta Prior Parameters for Bayesian estimation
#' @details  Evaluation of Bayesian Highest Probability Density (HPD) and two tailed
#' intervals using coverage probability, root mean square statistic, and the proportion
#' of proportion lies within the desired level of coverage for the Beta - Binomial
#' conjugate prior model for the probability of success \code{p}.
#' @return A dataframe with
#'  \item{method}{  Both Quantile and HPD method results are returned}
#'  \item{MeanCP}{  Coverage Probability}
#'  \item{MinCP }{  Minimum coverage probability}
#'  \item{RMSE_N }{   Root Mean Square Error from nominal size}
#'  \item{RMSE_M }{   Root Mean Square Error for Coverage Probability}
#'  \item{RMSE_MI }{  Root Mean Square Error for minimum coverage probability}
#'  \item{tol }{   Required tolerance for coverage probability}
#' @family Basic coverage probability methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; a=1;b=1; t1=0.93;t2=0.97;a1=1;a2=1
#' covpBA(n,alp,a,b,t1,t2,a1,a2)
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
#8.BAYESIAN
covpBA<-function(n,alp,a,b,t1,t2,a1,a2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if (missing(a1)) stop("'a1' is missing")
  if (missing(a2)) stop("'a2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")
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
cpBAQ=matrix(0,k,s)
ctBAQ=matrix(0,k,s)							#Cover Pbty quantity in sum
cppBAQ=0								#Coverage probabilty
RMSE_N1=0
RMSE_M1=0
RMSE_Mi1=0
ctr=0

cpBAH=matrix(0,k,s)
ctBAH=matrix(0,k,s)							#Cover Pbty quantity in sum
cppBAH=0								#Coverage probabilty
RMSE_N1H=0
RMSE_M1H=0
RMSE_Mi1H=0
ctrH=0

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
####COVERAGE PROBABILITIES
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
if(hp[j] > LBAQ[i] && hp[j] < UBAQ[i])
{
cpBAQ[i,j]=dbinom(i-1, n,hp[j])
ctBAQ[i,j]=1
}
if(hp[j] > LBAH[i] && hp[j] < UBAH[i])
{
cpBAH[i,j]=dbinom(i-1, n,hp[j])
ctBAH[i,j]=1
}

}
cppBAQ[j]=sum(cpBAQ[,j])						#Coverage Probability
RMSE_N1[j]=(cppBAQ[j]-(1-alp))^2			#Root mean Square from nominal size
if(t1<cppBAQ[j]&&cppBAQ[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined

cppBAH[j]=sum(cpBAH[,j])						#Coverage Probability
RMSE_N1H[j]=(cppBAH[j]-(1-alp))^2			#Root mean Square from nominal size
if(t1<cppBAH[j]&&cppBAH[j]<t2) ctrH=ctrH+1		#tolerance for cov prob - user defined

}
#CPBAQ=data.frame(hp,cppBAQ)
mcpBAQ=mean(cppBAQ)
micpBAQ=min(cppBAQ)					#Mean Cov Prob
RMSE_N=sqrt(mean(RMSE_N1))

#CPBAH=data.frame(hp,cppBAH)
mcpBAH=mean(cppBAH)
micpBAH=min(cppBAH)					#Mean Cov Prob
RMSE_NH=sqrt(mean(RMSE_N1H))


#Root mean Square from min and mean CoPr
for (j in 1:s)
{
RMSE_M1[j]=(cppBAQ[j]-mcpBAQ)^2
RMSE_Mi1[j]=(cppBAQ[j]-micpBAQ)^2
}
RMSE_M=sqrt(mean(RMSE_M1))
RMSE_MI=sqrt(mean(RMSE_Mi1))
tol=100*ctr/s

for (j in 1:s)
{
RMSE_M1H[j]=(cppBAH[j]-mcpBAQ)^2
RMSE_Mi1H[j]=(cppBAH[j]-micpBAQ)^2
}
RMSE_MH=sqrt(mean(RMSE_M1H))
RMSE_MIH=sqrt(mean(RMSE_Mi1H))
tolH=ctrH/s


TB=data.frame(Method="Quantile",MeanCP=mcpBAQ,MinCP=micpBAQ,RMSE_N=RMSE_N,RMSE_M=RMSE_M,RMSE_MI=RMSE_MI,tol=tol)
HPD=data.frame(Method="HPD",MeanCP=mcpBAH,MinCP=micpBAH,RMSE_N=RMSE_NH,RMSE_M=RMSE_MH,RMSE_MI=RMSE_MIH,tol=tolH)
ndf=rbind(TB,HPD)
return(ndf)
}

#############################################################################################################
#' Coverage Probability using 6 base methods
#' (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  The  Coverage Probability of 6 base methods
#' (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)
#' @return A dataframe with
#'  \item{method}{  method name}
#'  \item{MeanCP}{  Coverage Probability}
#'  \item{MinCP }{  Minimum coverage probability}
#'  \item{RMSE_N }{   Root Mean Square Error from nominal size}
#'  \item{RMSE_M }{   Root Mean Square Error for Coverage Probability}
#'  \item{RMSE_MI }{  Root Mean Square Error for minimum coverage probability}
#'  \item{tol }{   Required tolerance for coverage probability}
#' @family Basic coverage probability methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; a=1;b=1; t1=0.93;t2=0.97
#' covpAll(n,alp,a,b,t1,t2)
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
##### 9.All methods - Coverage Probability
covpAll<-function(n,alp,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")

  #### Calling functions and creating df
  Waldcovp.df    = covpWD(n,alp,a,b,t1,t2)
  ArcSinecovp.df = covpAS(n,alp,a,b,t1,t2)
  LRcovp.df      = covpLR(n,alp,a,b,t1,t2)
  Scorecovp.df   = covpSC(n,alp,a,b,t1,t2)
  WaldLcovp.df   = covpLT(n,alp,a,b,t1,t2)
  AdWaldcovp.df  = covpTW(n,alp,a,b,t1,t2)


  Waldcovp.df$method    = as.factor("Wald")
  ArcSinecovp.df$method = as.factor("ArcSine")
  LRcovp.df$method      = as.factor("Lilelihood")
  WaldLcovp.df$method    = as.factor("WaldLogit")
  Scorecovp.df$method   = as.factor("Score")
  AdWaldcovp.df$method  = as.factor("Wald-T")

  Generic.1 = data.frame(method = Waldcovp.df$method, MeanCP=Waldcovp.df$mcpW, MinCP= Waldcovp.df$micpW, RMSE_N=Waldcovp.df$RMSE_N,RMSE_M=Waldcovp.df$RMSE_M,RMSE_MI=Waldcovp.df$RMSE_MI,tol=Waldcovp.df$tol)
  Generic.2 = data.frame(method = ArcSinecovp.df$method, MeanCP=ArcSinecovp.df$mcpA, MinCP= ArcSinecovp.df$micpA, RMSE_N=ArcSinecovp.df$RMSE_N,RMSE_M=ArcSinecovp.df$RMSE_M,RMSE_MI=ArcSinecovp.df$RMSE_MI,tol=ArcSinecovp.df$tol)
  Generic.3 = data.frame(method = LRcovp.df$method, MeanCP=LRcovp.df$mcpL, MinCP= LRcovp.df$micpL, RMSE_N=LRcovp.df$RMSE_N,RMSE_M=LRcovp.df$RMSE_M,RMSE_MI=LRcovp.df$RMSE_MI,tol=LRcovp.df$tol)
  Generic.4 = data.frame(method = Scorecovp.df$method, MeanCP=Scorecovp.df$mcpS, MinCP= Scorecovp.df$micpS, RMSE_N=Scorecovp.df$RMSE_N,RMSE_M=Scorecovp.df$RMSE_M,RMSE_MI=Scorecovp.df$RMSE_MI,tol=Scorecovp.df$tol)
  Generic.5 = data.frame(method = WaldLcovp.df$method, MeanCP=WaldLcovp.df$mcpLT, MinCP= WaldLcovp.df$micpLT, RMSE_N=WaldLcovp.df$RMSE_N,RMSE_M=WaldLcovp.df$RMSE_M,RMSE_MI=WaldLcovp.df$RMSE_MI,tol=WaldLcovp.df$tol)
  Generic.6 = data.frame(method = AdWaldcovp.df$method, MeanCP=AdWaldcovp.df$mcpTW, MinCP= AdWaldcovp.df$micpTW, RMSE_N=AdWaldcovp.df$RMSE_N,RMSE_M=AdWaldcovp.df$RMSE_M,RMSE_MI=AdWaldcovp.df$RMSE_MI,tol=AdWaldcovp.df$tol)

  Final.df= rbind(Generic.1,Generic.2,Generic.3,Generic.4,Generic.5, Generic.6)

  return(Final.df)
}
