#' Coverage Probability of Continuity corrected Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continiuty correction
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  Evaluation of Wald-type interval with continuity correction using coverage
#' probability, root mean square statistic, and the proportion of proportion lies
#' within the desired level of coverage
#' @return A dataframe with
#'  \item{mcpCW}{  Continuity corrected Wald Coverage Probability}
#'  \item{micpCW }{   Continuity corrected Wald minimum coverage probability}
#'  \item{RMSE_N }{   Root Mean Square Error from nominal size}
#'  \item{RMSE_M }{   Root Mean Square Error for Coverage Probability}
#'  \item{RMSE_MI }{  Root Mean Square Error for minimum coverage probability}
#'  \item{tol }{   Required tolerance for coverage probability}
#' @family Coverage probability for continuity corrected methods
#' @examples
#' n= 10; alp=0.05; c=1/(2*n);a=1;b=1; t1=0.93;t2=0.97
#' covpCWD(n,alp,c,a,b,t1,t2)
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
##### 1.CC-WALD-Coverage Probability
covpCWD<-function(n,alp,c,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(c) != "integer") & (class(c) != "numeric") || length(c) >1 || c<0 ) stop("'c' has to be positive")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")

####INPUT n
x=0:n
k=n+1

####INITIALIZATIONS
pCW=0
qCW=0
seCW=0
LCW=0
UCW=0
s=5000								#Simulation run to generate hypothetical p
cpCW=matrix(0,k,s)
ctCW=matrix(0,k,s)							#Cover Pbty quantity in sum
cppCW=0								#Coverage probabilty
RMSE_N1=0
RMSE_M1=0
RMSE_Mi1=0
ctr=0
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
####COVERAGE PROBABILITIES
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
if(hp[j] > LCW[i] && hp[j] < UCW[i])
{
cpCW[i,j]=dbinom(i-1, n,hp[j])
ctCW[i,j]=1
}
}
cppCW[j]=sum(cpCW[,j])
RMSE_N1[j]=(cppCW[j]-(1-alp))^2			#Root mean Square from nominal size
if(t1<cppCW[j]&&cppCW[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined
}
#CPCW=data.frame(hp,cppCW)
mcpCW=mean(cppCW)
micpCW=min(cppCW)					#Mean Cov Prob
RMSE_N=sqrt(mean(RMSE_N1))

#Root mean Square from min and mean CoPr
for (j in 1:s)
{
RMSE_M1[j]=(cppCW[j]-mcpCW)^2
RMSE_Mi1[j]=(cppCW[j]-micpCW)^2
}
RMSE_M=sqrt(mean(RMSE_M1))
RMSE_MI=sqrt(mean(RMSE_Mi1))
tol=100*ctr/s
return(data.frame(mcpCW,micpCW,RMSE_N,RMSE_M,RMSE_MI,tol))
}

#############################################################################################################
#' Coverage Probability of Continuity corrected Score method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continiuty correction
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  Evaluation of continuity corrected score test approach using coverage probability,
#'  root mean square statistic, and the proportion of proportion lies within the
#'  desired level of coverage
#' @return A dataframe with
#'  \item{mcpAS}{  Continuity corrected Score Coverage Probability}
#'  \item{micpAS }{   Continuity corrected Score minimum coverage probability}
#'  \item{RMSE_N }{   Root Mean Square Error from nominal size}
#'  \item{RMSE_M }{   Root Mean Square Error for Coverage Probability}
#'  \item{RMSE_MI }{  Root Mean Square Error for minimum coverage probability}
#'  \item{tol }{   Required tolerance for coverage probability}
#' @family Coverage probability for continuity corrected methods
#' @examples
#' n= 10; alp=0.05; c=1/(2*n); a=1;b=1; t1=0.93;t2=0.97
#' covpCSC(n,alp,c,a,b,t1,t2)
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
##### 2.CC-SCORE - Coverage Probability
covpCSC<-function(n,alp,c,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if (c<=0 || c>(1/(2*n))) stop("'c' has to be positive and less than or equal to 1/(2*n)")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")


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
s=5000								#Simulation run to generate hypothetical p
cpCS=matrix(0,k,s)
ctCS=matrix(0,k,s)							#Cover Pbty quantity in sum
cppCS=0								#Coverage probabilty
RMSE_N1=0
RMSE_M1=0
RMSE_Mi1=0
ctr=0

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
####COVERAGE PROBABILITIES
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
if(hp[j] > LCS[i] && hp[j] < UCS[i])
{
cpCS[i,j]=dbinom(i-1, n,hp[j])
ctCS[i,j]=1
}
}
cppCS[j]=sum(cpCS[,j])						#Coverage Probability
RMSE_N1[j]=(cppCS[j]-(1-alp))^2			#Root mean Square from nominal size
if(t1<cppCS[j]&&cppCS[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined
}
#CPCS=data.frame(hp,cppCS)
mcpCS=mean(cppCS)							#Mean Cov Prob
micpCS=min(cppCS)							#Mean Cov Prob
RMSE_N=sqrt(mean(RMSE_N1))

#Root mean Square from min and mean CoPr
for (j in 1:s)
{
RMSE_M1[j]=(cppCS[j]-mcpCS)^2
RMSE_Mi1[j]=(cppCS[j]-micpCS)^2
}
RMSE_M=sqrt(mean(RMSE_M1))
RMSE_MI=sqrt(mean(RMSE_Mi1))
tol=100*ctr/s
return(data.frame(mcpCS,micpCS,RMSE_N,RMSE_M,RMSE_MI,tol))
}


#############################################################################################################
#' Coverage Probability of Continuity corrected ArcSine method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continiuty correction
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  Evaluation of continuity corrected Wald-type interval for the arcsine
#' transformation of the parameter p using coverage probability, root mean square
#' statistic, and the proportion of proportion lies within the desired level of coverage
#' @return A dataframe with
#'  \item{mcpCA}{  Continuity corrected ArcSine Coverage Probability}
#'  \item{micpCA }{   Continuity corrected ArcSine minimum coverage probability}
#'  \item{RMSE_N }{   Root Mean Square Error from nominal size}
#'  \item{RMSE_M }{   Root Mean Square Error for Coverage Probability}
#'  \item{RMSE_MI }{  Root Mean Square Error for minimum coverage probability}
#'  \item{tol }{   Required tolerance for coverage probability}
#' @family Coverage probability for continuity corrected methods
#' @examples
#' n= 10; alp=0.05; c=1/(2*n); a=1;b=1; t1=0.93;t2=0.97
#' covpCAS(n,alp,c,a,b,t1,t2)
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
##### 3.CC-ARC SINE - Coverage Probability
covpCAS<-function(n,alp,c,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(c) != "integer") & (class(c) != "numeric") || length(c) >1 || c<0 ) stop("'c' has to be positive")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")

####INPUT n
x=0:n
k=n+1
####INITIALIZATIONS
pCA=0
qCA=0
seCA=0
LCA=0
UCA=0
s=5000								#Simulation run to generate hypothetical p
cpCA=matrix(0,k,s)
ctCA=matrix(0,k,s)							#Cover Pbty quantity in sum
cppCA=0									#Coverage probabilty
RMSE_N1=0
RMSE_M1=0
RMSE_Mi1=0
ctr=0

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
####COVERAGE PROBABILITIES
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
if(hp[j] > LCA[i] && hp[j] < UCA[i])
{
cpCA[i,j]=dbinom(i-1, n,hp[j])
ctCA[i,j]=1
}
}
cppCA[j]=sum(cpCA[,j])						#Coverage Probability
RMSE_N1[j]=(cppCA[j]-(1-alp))^2			#Root mean Square from nominal size
if(t1<cppCA[j]&&cppCA[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined
}
#CPCA=data.frame(hp,cppCA)
mcpCA=mean(cppCA)							#Mean Cov Prob
micpCA=min(cppCA)
RMSE_N=sqrt(mean(RMSE_N1))

#Root mean Square from min and mean CoPr
for (j in 1:s)
{
RMSE_M1[j]=(cppCA[j]-mcpCA)^2
RMSE_Mi1[j]=(cppCA[j]-micpCA)^2
}
RMSE_M=sqrt(mean(RMSE_M1))
RMSE_MI=sqrt(mean(RMSE_Mi1))
tol=100*ctr/s
return(data.frame(mcpCA,micpCA,RMSE_N,RMSE_M,RMSE_MI,tol))
}
#############################################################################################################
#' Coverage Probability of Continuity corrected Logit Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continiuty correction
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  Evaluation of continuity corrected Wald-type interval based on the logit
#' transformation of p using coverage probability, root mean square statistic,
#' and the proportion of proportion lies within the desired level of coverage
#' @return A dataframe with
#'  \item{mcpALT}{  Continuity corrected Logit Wald Coverage Probability}
#'  \item{micpALT }{   Continuity corrected Logit Wald minimum coverage probability}
#'  \item{RMSE_N }{   Root Mean Square Error from nominal size}
#'  \item{RMSE_M }{   Root Mean Square Error for Coverage Probability}
#'  \item{RMSE_MI }{  Root Mean Square Error for minimum coverage probability}
#'  \item{tol }{   Required tolerance for coverage probability}
#' @family Coverage probability for continuity corrected methods
#' @examples
#' n= 10; alp=0.05;c=1/(2*n); a=1;b=1; t1=0.93;t2=0.97
#' covpCLT(n,alp,c,a,b,t1,t2)
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
##### 4.LOGIT-WALD - Coverage Probability
covpCLT<-function(n,alp,c,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(c) != "integer") & (class(c) != "numeric") || length(c) >1 || c<0 ) stop("'c' has to be positive")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")

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
s=5000								#Simulation run to generate hypothetical p
cpCLT=matrix(0,k,s)
ctCLT=matrix(0,k,s)							#Cover Pbty quantity in sum
cppCLT=0								#Coverage probabilty
RMSE_N1=0
RMSE_M1=0
RMSE_Mi1=0
ctr=0

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
}
####COVERAGE PROBABILITIES
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
if(hp[j] > LCLT[i] && hp[j] < UCLT[i])
{
cpCLT[i,j]=dbinom(i-1, n,hp[j])
ctCLT[i,j]=1
}
}
cppCLT[j]=sum(cpCLT[,j])				#Coverage Probability
RMSE_N1[j]=(cppCLT[j]-(1-alp))^2			#Root mean Square from nominal size
if(t1<cppCLT[j]&&cppCLT[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined

}
#CPLT=data.frame(hp,cppCLT)
mcpCLT=mean(cppCLT)
micpCLT=min(cppCLT)					#Mean Cov Prob
RMSE_N=sqrt(mean(RMSE_N1))

#Root mean Square from min and mean CoPr
for (j in 1:s)
{
RMSE_M1[j]=(cppCLT[j]-mcpCLT)^2
RMSE_Mi1[j]=(cppCLT[j]-micpCLT)^2
}
RMSE_M=sqrt(mean(RMSE_M1))
RMSE_MI=sqrt(mean(RMSE_Mi1))
tol=100*ctr/s
return(data.frame(mcpCLT,micpCLT,RMSE_N,RMSE_M,RMSE_MI,tol))
}
#############################################################################################################
#' Coverage Probability of Continuity corrected Wald-T method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continiuty correction
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  Evaluation of approximate and continuity corrected method based on a
#' t_approximation of the standardized point estimator using coverage probability,
#' root mean square statistic, and the proportion of proportion lies within the
#' desired level of coverage
#' @return A dataframe with
#'  \item{mcpATW}{  Continuity corrected Wald-T Coverage Probability}
#'  \item{micpATW }{   Continuity corrected Wald-T minimum coverage probability}
#'  \item{RMSE_N }{   Root Mean Square Error from nominal size}
#'  \item{RMSE_M }{   Root Mean Square Error for Coverage Probability}
#'  \item{RMSE_MI }{  Root Mean Square Error for minimum coverage probability}
#'  \item{tol }{   Required tolerance for coverage probability}
#' @family Coverage probability for continuity corrected methods
#' @examples
#' n= 10; alp=0.05; c=1/(2*n);a=1;b=1; t1=0.93;t2=0.97
#' covpCTW(n,alp,c,a,b,t1,t2)
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
##### 5. WALD_t - Coverage Probability
covpCTW<-function(n,alp,c,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if ((class(c) != "integer") & (class(c) != "numeric") || length(c) >1 || c<0 ) stop("'c' has to be positive")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")

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
s=5000								#Simulation run to generate hypothetical p
cpCTW=matrix(0,k,s)
ctCTW=matrix(0,k,s)							#Cover Pbty quantity in sum
cppCTW=0
RMSE_N1=0
RMSE_M1=0
RMSE_Mi1=0
ctr=0

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
####COVERAGE PROBABILITIES
hp=sort(rbeta(s,a,b),decreasing = FALSE)	#HYPOTHETICAL "p"
for (j in 1:s)
{
for(i in 1:k)
{
if(hp[j] > LCTW[i] && hp[j] < UCTW[i])
{
cpCTW[i,j]=dbinom(i-1, n,hp[j])
ctCTW[i,j]=1
}
}
cppCTW[j]=sum(cpCTW[,j])						#Coverage Probability
RMSE_N1[j]=(cppCTW[j]-(1-alp))^2			#Root mean Square from nominal size
if(t1<cppCTW[j]&&cppCTW[j]<t2) ctr=ctr+1		#tolerance for cov prob - user defined
}
#CPCTW=data.frame(hp,cppCTW)
mcpCTW=mean(cppCTW)
micpCTW=min(cppCTW)					#Mean Cov Prob
RMSE_N=sqrt(mean(RMSE_N1))

#Root mean Square from min and mean CoPr
for (j in 1:s)
{
RMSE_M1[j]=(cppCTW[j]-mcpCTW)^2
RMSE_Mi1[j]=(cppCTW[j]-micpCTW)^2
}
RMSE_M=sqrt(mean(RMSE_M1))
RMSE_MI=sqrt(mean(RMSE_Mi1))
tol=100*ctr/s
return(data.frame(mcpCTW,micpCTW,RMSE_N,RMSE_M,RMSE_MI,tol))
}

#############################################################################################################
#' Coverage Probability for 5 continuity corrected methods (Wald, Wald-T, Score, Logit-Wald, ArcSine)
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param c - Continiuty correction
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @param t1 - Lower tolerance limit to check the spread of coverage Probability
#' @param t2 - Upper tolerance limit to check the spread of coverage Probability
#' @details  The  Coverage Probability of 5 continuity corrected methods (Wald, Wald-T, Score, Logit-Wald, ArcSine) for \code{n} given \code{alp}, \code{h}, \code{a}, \code{b}, \code{t1} and  \code{t2} using all the methods
#' @return A dataframe with
#'  \item{method}{  Method name}
#'  \item{MeanCP}{  Coverage Probability}
#'  \item{MinCP }{  Minimum coverage probability}
#'  \item{RMSE_N }{   Root Mean Square Error from nominal size}
#'  \item{RMSE_M }{   Root Mean Square Error for Coverage Probability}
#'  \item{RMSE_MI }{  Root Mean Square Error for minimum coverage probability}
#'  \item{tol }{   Required tolerance for coverage probability}
#' @family Coverage probability for continuity corrected methods
#' @examples
#' \dontrun{
#' n= 10; alp=0.05; c=1/(2*n);a=1;b=1; t1=0.93;t2=0.97
#' covpCAll(n,alp,c,a,b,t1,t2)
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
##### 9.All methods - Coverage Probability=5 continuity corrected methods (Wald, Wald-T, Score, Logit-Wald, ArcSine)
covpCAll<-function(n,alp,c,a,b,t1,t2)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(c)) stop("'c' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(t1)) stop("'t1' is missing")
  if (missing(t2)) stop("'t2' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp) >1) stop("'alpha' has to be between 0 and 1")
  if (c<=0 || c>(1/(2*n))) stop("'c' has to be positive and less than or equal to 1/(2*n)")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a)>1 || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b)>1 || b<0  ) stop("'b' has to be greater than or equal to 0")
  if (t1>t2) stop(" t1 has to be lesser than t2")
  if ((class(t1) != "integer") & (class(t1) != "numeric") || length(t1)>1 || t1<0 || t1>1 ) stop("'t1' has to be between 0 and 1")
  if ((class(t2) != "integer") & (class(t2) != "numeric") || length(t2)>1 || t2<0 || t2>1 ) stop("'t2' has to be between 0 and 1")

  #### Calling functions and creating df
  WaldcovpA.df    = covpCWD(n,alp,c,a,b,t1,t2)
  ArcSinecovpA.df = covpCAS(n,alp,c,a,b,t1,t2)
  ScorecovpA.df   = covpCSC(n,alp,c,a,b,t1,t2)
  WaldLcovpA.df   = covpCLT(n,alp,c,a,b,t1,t2)
  AdWaldcovpA.df  = covpCTW(n,alp,c,a,b,t1,t2)


  WaldcovpA.df$method    = as.factor("CC-Wald")
  ArcSinecovpA.df$method = as.factor("CC-ArcSine")
  WaldLcovpA.df$method    = as.factor("CC-Logit-Wald")
  ScorecovpA.df$method   = as.factor("CC-Score")
  AdWaldcovpA.df$method  = as.factor("CC-Wald-T")

  Generic.1 = data.frame(method = WaldcovpA.df$method, MeanCP=WaldcovpA.df$mcpCW, MinCP= WaldcovpA.df$micpCW, RMSE_N=WaldcovpA.df$RMSE_N,RMSE_M=WaldcovpA.df$RMSE_M,RMSE_MI=WaldcovpA.df$RMSE_MI,tol=WaldcovpA.df$tol)
  Generic.2 = data.frame(method = ArcSinecovpA.df$method, MeanCP=ArcSinecovpA.df$mcpCA, MinCP= ArcSinecovpA.df$micpCA, RMSE_N=ArcSinecovpA.df$RMSE_N,RMSE_M=ArcSinecovpA.df$RMSE_M,RMSE_MI=ArcSinecovpA.df$RMSE_MI,tol=ArcSinecovpA.df$tol)
  Generic.4 = data.frame(method = ScorecovpA.df$method, MeanCP=ScorecovpA.df$mcpCS, MinCP= ScorecovpA.df$micpCS, RMSE_N=ScorecovpA.df$RMSE_N,RMSE_M=ScorecovpA.df$RMSE_M,RMSE_MI=ScorecovpA.df$RMSE_MI,tol=ScorecovpA.df$tol)
  Generic.5 = data.frame(method = WaldLcovpA.df$method, MeanCP=WaldLcovpA.df$mcpCLT, MinCP= WaldLcovpA.df$micpCLT, RMSE_N=WaldLcovpA.df$RMSE_N,RMSE_M=WaldLcovpA.df$RMSE_M,RMSE_MI=WaldLcovpA.df$RMSE_MI,tol=WaldLcovpA.df$tol)
  Generic.6 = data.frame(method = AdWaldcovpA.df$method, MeanCP=AdWaldcovpA.df$mcpCTW, MinCP= AdWaldcovpA.df$micpCTW, RMSE_N=AdWaldcovpA.df$RMSE_N,RMSE_M=AdWaldcovpA.df$RMSE_M,RMSE_MI=AdWaldcovpA.df$RMSE_MI,tol=AdWaldcovpA.df$tol)

  Final.df= rbind(Generic.1,Generic.2,Generic.4,Generic.5, Generic.6)

  return(Final.df)
}
