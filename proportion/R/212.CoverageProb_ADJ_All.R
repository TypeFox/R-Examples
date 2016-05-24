#' Coverage Probability of Adjusted Wald method for given n
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  Evaluation of adjusted Wald-type interval using coverage probability,
#' root mean square statistic, and the proportion of proportion lies within the
#' desired level of coverage
#' @return A dataframe with
#'  \item{mcpAW}{  Adjusted Wald Coverage Probability}
#'  \item{micpAW }{   Adjusted Wald minimum coverage probability}
#'  \item{RMSE_N }{   Root Mean Square Error from nominal size}
#'  \item{RMSE_M }{   Root Mean Square Error for Coverage Probability}
#'  \item{RMSE_MI }{  Root Mean Square Error for minimum coverage probability}
#'  \item{tol }{   Required tolerance for coverage probability}
#' @family Coverage probability of adjusted methods
#' @examples
#' n= 10; alp=0.05; h=2;a=1;b=1; t1=0.93;t2=0.97
#' covpAWD(n,alp,h,a,b,t1,t2)
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
##### 1.ADJUSTED WALD-Coverage Probability
covpAWD<-function(n,alp,h,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h)>1 || h<0  ) stop("'h' has to be greater than or equal to 0")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")

####INPUT n
x=0:n
k=n+1
y=x+h
n1=n+(2*h)
####INITIALIZATIONS
pAW=0
qAW=0
seAW=0
LAW=0
UAW=0
s=5000								#Simulation run to generate hypothetical p
cpAW=matrix(0,k,s)
ctAW=matrix(0,k,s)							#Cover Pbty quantity in sum
cppAW=0								#Coverage probabilty
RMSE_N1=0
RMSE_M1=0
RMSE_Mi1=0
ctr=0
###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#WALD METHOD
for(i in 1:k)
{
pAW[i]=y[i]/n1
qAW[i]=1-pAW[i]
seAW[i]=sqrt(pAW[i]*qAW[i]/n1)
LAW[i]=pAW[i]-(cv*seAW[i])
UAW[i]=pAW[i]+(cv*seAW[i])
if(LAW[i]<0) LAW[i]=0
if(UAW[i]>1) UAW[i]=1
}
####COVERAGE PROBABILITIES
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
if(hp[j] > LAW[i] && hp[j] < UAW[i])
{
cpAW[i,j]=dbinom(i-1, n,hp[j])
ctAW[i,j]=1
}
}
cppAW[j]=sum(cpAW[,j])
RMSE_N1[j]=(cppAW[j]-(1-alp))^2			#Root mean Square from nominal size
if(t1<cppAW[j]&&cppAW[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined
}
#CPAW=data.frame(hp,cppAW)
mcpAW=mean(cppAW)
micpAW=min(cppAW)					#Mean Cov Prob
RMSE_N=sqrt(mean(RMSE_N1))

#Root mean Square from min and mean CoPr
for (j in 1:s)
{
RMSE_M1[j]=(cppAW[j]-mcpAW)^2
RMSE_Mi1[j]=(cppAW[j]-micpAW)^2
}
RMSE_M=sqrt(mean(RMSE_M1))
RMSE_MI=sqrt(mean(RMSE_Mi1))
tol=100*ctr/s
return(data.frame(mcpAW,micpAW,RMSE_N,RMSE_M,RMSE_MI,tol))
}

#############################################################################################################
#' Coverage Probability of Adjusted Score method  for given n
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  Evaluation of adjusted score test approach using coverage probability,
#' root mean square statistic, and the proportion of proportion lies within the
#' desired level of coverage
#' @return A dataframe with
#'  \item{mcpAS}{  Adjusted Score Coverage Probability}
#'  \item{micpAS }{   Adjusted Score minimum coverage probability}
#'  \item{RMSE_N }{   Root Mean Square Error from nominal size}
#'  \item{RMSE_M }{   Root Mean Square Error for Coverage Probability}
#'  \item{RMSE_MI }{  Root Mean Square Error for minimum coverage probability}
#'  \item{tol }{   Required tolerance for coverage probability}
#' @family Coverage probability of adjusted methods
#' @examples
#' n= 10; alp=0.05; h=2; a=1;b=1; t1=0.93;t2=0.97
#' covpASC(n,alp,h,a,b,t1,t2)
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
##### 2.ADJUSTED SCORE - Coverage Probability
covpASC<-function(n,alp,h,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h)>1 || h<0  ) stop("'h' has to be greater than or equal to 0")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")

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
s=1000								#Simulation run to generate hypothetical p
cpAS=matrix(0,k,s)
ctAS=matrix(0,k,s)							#Cover Pbty quantity in sum
cppAS=0							#Coverage probabilty
RMSE_N1=0
RMSE_M1=0
RMSE_Mi1=0
ctr=0

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
####COVERAGE PROBABILITIES
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
if(hp[j] > LAS[i] && hp[j] < UAS[i])
{
cpAS[i,j]=dbinom(i-1, n,hp[j])
ctAS[i,j]=1
}
}
cppAS[j]=sum(cpAS[,j])						#Coverage Probability
RMSE_N1[j]=(cppAS[j]-(1-alp))^2			#Root mean Square from nominal size
if(t1<cppAS[j]&&cppAS[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined
}
#CPAS=data.frame(hp,cppAS)
mcpAS=mean(cppAS)							#Mean Cov Prob
micpAS=min(cppAS)							#Mean Cov Prob
RMSE_N=sqrt(mean(RMSE_N1))

#Root mean Square from min and mean CoPr
for (j in 1:s)
{
RMSE_M1[j]=(cppAS[j]-mcpAS)^2
RMSE_Mi1[j]=(cppAS[j]-micpAS)^2
}
RMSE_M=sqrt(mean(RMSE_M1))
RMSE_MI=sqrt(mean(RMSE_Mi1))
tol=100*ctr/s
return(data.frame(mcpAS,micpAS,RMSE_N,RMSE_M,RMSE_MI,tol))
}


#############################################################################################################
#' Coverage Probability of Adjusted ArcSine method  for given n
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  Evaluation of adjusted Wald-type interval for the arcsine transformation of
#' the parameter p using coverage probability, root mean square statistic, and the
#' proportion of proportion lies within the desired level of coverage
#' @return A dataframe with
#'  \item{mcpAA}{  Adjusted ArcSine Coverage Probability}
#'  \item{micpAA }{   Adjusted ArcSine minimum coverage probability}
#'  \item{RMSE_N }{   Root Mean Square Error from nominal size}
#'  \item{RMSE_M }{   Root Mean Square Error for Coverage Probability}
#'  \item{RMSE_MI }{  Root Mean Square Error for minimum coverage probability}
#'  \item{tol }{   Required tolerance for coverage probability}
#' @family Coverage probability of adjusted methods
#' @examples
#' n= 10; alp=0.05; h=2; a=1;b=1; t1=0.93;t2=0.97
#' covpAAS(n,alp,h,a,b,t1,t2)
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
##### 3.ADJUSTED ARC SINE - Coverage Probability
covpAAS<-function(n,alp,h,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h)>1 || h<0  ) stop("'h' has to be greater than or equal to 0")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")

####INPUT n
x=0:n
k=n+1
y=x+h
n1=n+(2*h)
####INITIALIZATIONS
pAA=0
qAA=0
seAA=0
LAA=0
UAA=0
s=5000								#Simulation run to generate hypothetical p
cpAA=matrix(0,k,s)
ctAA=matrix(0,k,s)							#Cover Pbty quantity in sum
cppAA=0								#Coverage probabilty
RMSE_N1=0
RMSE_M1=0
RMSE_Mi1=0
ctr=0

###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#ADJUSTED ARC-SINE METHOD
for(i in 1:k)
{
pAA[i]=y[i]/n1
qAA[i]=1-pAA[i]
seAA[i]=cv/sqrt(4*n1)
LAA[i]=(sin(asin(sqrt(pAA[i]))-seAA[i]))^2
UAA[i]=(sin(asin(sqrt(pAA[i]))+seAA[i]))^2
if(LAA[i]<0) LAA[i]=0
if(UAA[i]>1) UAA[i]=1
}
####COVERAGE PROBABILITIES
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
if(hp[j] > LAA[i] && hp[j] < UAA[i])
{
cpAA[i,j]=dbinom(i-1, n,hp[j])
ctAA[i,j]=1
}
}
cppAA[j]=sum(cpAA[,j])						#Coverage Probability
RMSE_N1[j]=(cppAA[j]-(1-alp))^2			#Root mean Square from nominal size
if(t1<cppAA[j]&&cppAA[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined
}
#CPAA=data.frame(hp,cppAA)
mcpAA=mean(cppAA)							#Mean Cov Prob
micpAA=min(cppAA)
RMSE_N=sqrt(mean(RMSE_N1))

#Root mean Square from min and mean CoPr
for (j in 1:s)
{
RMSE_M1[j]=(cppAA[j]-mcpAA)^2
RMSE_Mi1[j]=(cppAA[j]-micpAA)^2
}
RMSE_M=sqrt(mean(RMSE_M1))
RMSE_MI=sqrt(mean(RMSE_Mi1))
tol=100*ctr/s
return(data.frame(mcpAA,micpAA,RMSE_N,RMSE_M,RMSE_MI,tol))
}
#############################################################################################################
#' Coverage Probability of Adjusted Logit Wald method  for given n
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  Evaluation of adjusted Wald-type interval based on the
#' logit transformation of \code{p} using coverage probability, root mean square statistic,
#' and the proportion of proportion lies within the desired level of coverage.
#' @return A dataframe with
#'  \item{mcpALT}{  Adjusted Logit Wald Coverage Probability}
#'  \item{micpALT }{   Adjusted Logit Wald minimum coverage probability}
#'  \item{RMSE_N }{   Root Mean Square Error from nominal size}
#'  \item{RMSE_M }{   Root Mean Square Error for Coverage Probability}
#'  \item{RMSE_MI }{  Root Mean Square Error for minimum coverage probability}
#'  \item{tol }{   Required tolerance for coverage probability}
#' @family Coverage probability of adjusted methods
#' @examples
#' n= 10; alp=0.05;h=2; a=1;b=1; t1=0.93;t2=0.97
#' covpALT(n,alp,h,a,b,t1,t2)
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
##### 4.ADJUSTED LOGIT-WALD - Coverage Probability
covpALT<-function(n,alp,h,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h)>1 || h<0  ) stop("'h' has to be greater than or equal to 0")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")

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
s=5000								#Simulation run to generate hypothetical p
cpALT=matrix(0,k,s)
ctALT=matrix(0,k,s)							#Cover Pbty quantity in sum
cppALT=0								#Coverage probabilty
RMSE_N1=0
RMSE_M1=0
RMSE_Mi1=0
ctr=0

###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#ADJUSTED LOGIT-WALD METHOD
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
####COVERAGE PROBABILITIES
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
if(hp[j] > LALT[i] && hp[j] < UALT[i])
{
cpALT[i,j]=dbinom(i-1, n,hp[j])
ctALT[i,j]=1
}
}
cppALT[j]=sum(cpALT[,j])				#Coverage Probability
RMSE_N1[j]=(cppALT[j]-(1-alp))^2			#Root mean Square from nominal size
if(t1<cppALT[j]&&cppALT[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined

}
#CPALT=data.frame(hp,cppALT)
mcpALT=mean(cppALT)
micpALT=min(cppALT)					#Mean Cov Prob
RMSE_N=sqrt(mean(RMSE_N1))

#Root mean Square from min and mean CoPr
for (j in 1:s)
{
RMSE_M1[j]=(cppALT[j]-mcpALT)^2
RMSE_Mi1[j]=(cppALT[j]-micpALT)^2
}
RMSE_M=sqrt(mean(RMSE_M1))
RMSE_MI=sqrt(mean(RMSE_Mi1))
tol=100*ctr/s
return(data.frame(mcpALT,micpALT,RMSE_N,RMSE_M,RMSE_MI,tol))
}
#############################################################################################################
#' Coverage Probability of Adjusted Wald-T method  for given n
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  Evaluation of approximate and adjusted method based on a t_approximation of
#' the standardized point estimator using coverage probability, root mean square statistic,
#' and the proportion of proportion lies within the desired level of coverage
#' @return A dataframe with
#'  \item{mcpATW}{  Adjusted Wald-T Coverage Probability}
#'  \item{micpATW }{   Adjusted Wald-T minimum coverage probability}
#'  \item{RMSE_N }{   Root Mean Square Error from nominal size}
#'  \item{RMSE_M }{   Root Mean Square Error for Coverage Probability}
#'  \item{RMSE_MI }{  Root Mean Square Error for minimum coverage probability}
#'  \item{tol }{   Required tolerance for coverage probability}
#' @family Coverage probability of adjusted methods
#' @examples
#' n= 10; alp=0.05; h=2;a=1;b=1; t1=0.93;t2=0.97
#' covpATW(n,alp,h,a,b,t1,t2)
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
##### 5. ADJUSTED WALD_t - Coverage Probability
covpATW<-function(n,alp,h,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h)>1 || h<0  ) stop("'h' has to be greater than or equal to 0")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")

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
s=5000								#Simulation run to generate hypothetical p
cpATW=matrix(0,k,s)
ctATW=matrix(0,k,s)							#Cover Pbty quantity in sum
cppATW=0
RMSE_N1=0
RMSE_M1=0
RMSE_Mi1=0
ctr=0

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
####COVERAGE PROBABILITIES
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
if(hp[j] > LATW[i] && hp[j] < UATW[i])
{
cpATW[i,j]=dbinom(i-1, n,hp[j])
ctATW[i,j]=1
}
}
cppATW[j]=sum(cpATW[,j])						#Coverage Probability
RMSE_N1[j]=(cppATW[j]-(1-alp))^2			#Root mean Square from nominal size
if(t1<cppATW[j]&&cppATW[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined
}
#CPATW=data.frame(hp,cppATW)
mcpATW=mean(cppATW)
micpATW=min(cppATW)					#Mean Cov Prob
RMSE_N=sqrt(mean(RMSE_N1))

#Root mean Square from min and mean CoPr
for (j in 1:s)
{
RMSE_M1[j]=(cppATW[j]-mcpATW)^2
RMSE_Mi1[j]=(cppATW[j]-micpATW)^2
}
RMSE_M=sqrt(mean(RMSE_M1))
RMSE_MI=sqrt(mean(RMSE_Mi1))
tol=100*ctr/s
return(data.frame(mcpATW,micpATW,RMSE_N,RMSE_M,RMSE_MI,tol))
}

#############################################################################################################
#'    Coverage Probability of Adjusted Likelihood  method for given n
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  Evaluation of adjusted Likelihood ratio limits using coverage probability,
#' root mean square statistic, and the proportion of proportion lies within the desired
#' level of coverage
#' @return A dataframe with
#'  \item{mcpAL}{   Adjusted Likelihood  Coverage Probability}
#'  \item{micpAL }{    Adjusted Likelihood  minimum coverage probability}
#'  \item{RMSE_N }{   Root Mean Square Error from nominal size}
#'  \item{RMSE_M }{   Root Mean Square Error for Coverage Probability}
#'  \item{RMSE_MI }{  Root Mean Square Error for minimum coverage probability}
#'  \item{tol }{   Required tolerance for coverage probability}
#' @family Coverage probability of adjusted methods
#' @examples
#' n= 10; alp=0.05; h=2;a=1;b=1; t1=0.93;t2=0.97
#' covpALR(n,alp,h,a,b,t1,t2)
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
##### 6.ADJUSTED LIKELIHOOD RATIO - Coverage Probability
covpALR<-function(n,alp,h,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h) >1|| h<0  || !(h%%1 ==0)) stop("'h' has to be an integer greater than or equal to 0")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")

####INPUT n
y=0:n
k=n+1
y1=y+h
n1=n+(2*h)
####INITIALIZATIONS
mle=0
cutoff=0
LAL=0
UAL=0
s=5000								#Simulation run to generate hypothetical p
cpAL=matrix(0,k,s)
ctAL=matrix(0,k,s)							#Cover Pbty quantity in sum
cppAL=0								#Coverage probabilty
RMSE_N1=0
RMSE_M1=0
RMSE_Mi1=0
ctr=0

###CRITICAL VALUES
cv=qnorm(1-(alp/2), mean = 0, sd = 1)
#ADJUSTED LIKELIHOOD-RATIO METHOD
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
####COVERAGE PROBABILITIES
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
if(hp[j] > LAL[i] && hp[j] < UAL[i])
{
cpAL[i,j]=dbinom(i-1, n,hp[j])
ctAL[i,j]=1
}
}
cppAL[j]=sum(cpAL[,j])						#Coverage Probability
RMSE_N1[j]=(cppAL[j]-(1-alp))^2			#Root mean Square from nominal size
if(t1<cppAL[j]&&cppAL[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined
}
#CPAL=data.frame(hp,cppAL)
mcpAL=mean(cppAL)
micpAL=min(cppAL)					#Mean Cov Prob
RMSE_N=sqrt(mean(RMSE_N1))

#Root mean Square from min and mean CoPr
for (j in 1:s)
{
RMSE_M1[j]=(cppAL[j]-mcpAL)^2
RMSE_Mi1[j]=(cppAL[j]-micpAL)^2
}
RMSE_M=sqrt(mean(RMSE_M1))
RMSE_MI=sqrt(mean(RMSE_Mi1))
tol=100*ctr/s
return(data.frame(mcpAL,micpAL,RMSE_N,RMSE_M,RMSE_MI,tol))
}

#############################################################################################################
#' Coverage Probability for 6 adjusted methods
#' (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  Calculates the  Coverage Probability  for 6 adjusted methods
#' (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)
#' @return A dataframe with
#'  \item{method}{  Method name}
#'  \item{MeanCP}{  Coverage Probability}
#'  \item{MinCP }{  Minimum coverage probability}
#'  \item{RMSE_N }{   Root Mean Square Error from nominal size}
#'  \item{RMSE_M }{   Root Mean Square Error for Coverage Probability}
#'  \item{RMSE_MI }{  Root Mean Square Error for minimum coverage probability}
#'  \item{tol }{   Required tolerance for coverage probability}
#' @family Coverage probability of adjusted methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; h=2;a=1;b=1; t1=0.93;t2=0.97
#' covpAAll(n,alp,h,a,b,t1,t2)
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
##### 9.All methods - Coverage Probability
covpAAll<-function(n,alp,h,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h) >1|| h<0  || !(h%%1 ==0)) stop("'h' has to be an integer greater than or equal to 0")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")

  #### Calling functions and creating df
  Waldcovp.df    = covpAWD(n,alp,h,a,b,t1,t2)
  ArcSinecovp.df = covpAAS(n,alp,h,a,b,t1,t2)
  LRcovp.df      = covpALR(n,alp,h,a,b,t1,t2)
  Scorecovp.df   = covpASC(n,alp,h,a,b,t1,t2)
  WaldLcovp.df   = covpALT(n,alp,h,a,b,t1,t2)
  AdWaldcovp.df  = covpATW(n,alp,h,a,b,t1,t2)


  Waldcovp.df$method    = as.factor("Adj-Wald")
  ArcSinecovp.df$method = as.factor("Adj-ArcSine")
  LRcovp.df$method      = as.factor("Adj-Lilelihood")
  WaldLcovp.df$method    = as.factor("Adj-Logit-Wald")
  Scorecovp.df$method   = as.factor("Adj-Score")
  AdWaldcovp.df$method  = as.factor("Adj-Wald-T")

  Generic.1 = data.frame(method = Waldcovp.df$method, MeanCP=Waldcovp.df$mcpAW, MinCP= Waldcovp.df$micpAW, RMSE_N=Waldcovp.df$RMSE_N,RMSE_M=Waldcovp.df$RMSE_M,RMSE_MI=Waldcovp.df$RMSE_MI,tol=Waldcovp.df$tol)
  Generic.2 = data.frame(method = ArcSinecovp.df$method, MeanCP=ArcSinecovp.df$mcpAA, MinCP= ArcSinecovp.df$micpAA, RMSE_N=ArcSinecovp.df$RMSE_N,RMSE_M=ArcSinecovp.df$RMSE_M,RMSE_MI=ArcSinecovp.df$RMSE_MI,tol=ArcSinecovp.df$tol)
  Generic.3 = data.frame(method = LRcovp.df$method, MeanCP=LRcovp.df$mcpAL, MinCP= LRcovp.df$micpAL, RMSE_N=LRcovp.df$RMSE_N,RMSE_M=LRcovp.df$RMSE_M,RMSE_MI=LRcovp.df$RMSE_MI,tol=LRcovp.df$tol)
  Generic.4 = data.frame(method = Scorecovp.df$method, MeanCP=Scorecovp.df$mcpAS, MinCP= Scorecovp.df$micpAS, RMSE_N=Scorecovp.df$RMSE_N,RMSE_M=Scorecovp.df$RMSE_M,RMSE_MI=Scorecovp.df$RMSE_MI,tol=Scorecovp.df$tol)
  Generic.5 = data.frame(method = WaldLcovp.df$method, MeanCP=WaldLcovp.df$mcpALT, MinCP= WaldLcovp.df$micpALT, RMSE_N=WaldLcovp.df$RMSE_N,RMSE_M=WaldLcovp.df$RMSE_M,RMSE_MI=WaldLcovp.df$RMSE_MI,tol=WaldLcovp.df$tol)
  Generic.6 = data.frame(method = AdWaldcovp.df$method, MeanCP=AdWaldcovp.df$mcpATW, MinCP= AdWaldcovp.df$micpATW, RMSE_N=AdWaldcovp.df$RMSE_N,RMSE_M=AdWaldcovp.df$RMSE_M,RMSE_MI=AdWaldcovp.df$RMSE_MI,tol=AdWaldcovp.df$tol)

  Final.df= rbind(Generic.1,Generic.2,Generic.3,Generic.4,Generic.5, Generic.6)

  return(Final.df)
}
