#' Permutation Using the Maximum Cross Correlation Method
#' 
#' The permutation is determined by grouping the components of a multivariate series X into \eqn{q} groups, where \eqn{q} and the cardinal numbers of those groups are also unknown.
#'
#' See Chang et al. (2014) for the permutation step and more information. 
#'
#' @param X a data matrix used to find the grouping mechanism with \eqn{n} rows and \eqn{p} columns, where \eqn{n} is the sample size and \eqn{p} is the dimension of the time series. 
#' @param Vol logical. If \code{FALSE} (the default), then prewhiten each series by fitting a univariate AR model with
#'          the order between 0 and 5 determined by AIC. If \code{TRUE}, then prewhiten each volatility process using GARCH(1,1) model.                 
#' @param m a positive constant used to calculate the maximum cross correlation over the lags between \eqn{-m} and \eqn{m}. If \eqn{m} is not specified, the default constant \eqn{10*log10(n/p)} 
#'          will be used. 

#' @return An object of class "permutationMax" is a list containing the following components:
#'
#' \item{NoGroups}{number of groups with at least two components series}
#' \item{Nos_of_Members}{number of members in each of groups with at least two members} 		   
#' \item{Groups}{indices of components in each of groups with at least two members}
#' \item{maxcorr}{maximum correlation (over lags) of \eqn{p(p-1)/2} pairs in descending order}
#' \item{corrRatio}{ratios of successive values from maxcorr}
#' \item{NoConnectedPairs}{number of connected pairs}
#' \item{Xpre}{the prewhitened data with \eqn{n-R} rows and \eqn{p} columns}
#' @note This is the second step for segmentation by grouping the transformed time series. The first step is to seek for a contemporaneous linear transformation of the original series, see \code{\link{segmentTS}}.      
 
#' @author Jinyuan Chang, Bin Guo and Qiwei Yao
#' @references Chang, J., Guo, B. and Yao, Q. (2014). Segmenting Multiple Time Series by Contemporaneous Linear Transformation: PCA for Time Series. Available at \url{http://arxiv.org/abs/1410.2323}
#' @seealso \code{\link{segmentTS}}, \code{\link{permutationFDR}}
#' @examples
#' ## Example 1 (Example 5 of Chang et al.(2014)).
#' ## p=6, x_t consists of 3 independent subseries with 3, 2 and 1 components.    
#'
#' p=6;n=1500
#' # Generate x_t
#' X=mat.or.vec(p,n)
#' x=arima.sim(model=list(ar=c(0.5, 0.3), ma=c(-0.9, 0.3, 1.2,1.3)),n=n+2,sd=1)
#' for(i in 1:3) X[i,]=x[i:(n+i-1)]  
#' x=arima.sim(model=list(ar=c(0.8,-0.5),ma=c(1,0.8,1.8) ),n=n+1,sd=1)
#' for(i in 4:5) X[i,]=x[(i-3):(n+i-4)]   
#' x=arima.sim(model=list(ar=c(-0.7, -0.5), ma=c(-1, -0.8)),n=n,sd=1)
#' X[6,]=x
#' # Generate y_t 
#' A=matrix(runif(p*p, -3, 3), ncol=p)
#' Y=A%*%X  
#' Y=t(Y)
#' Trans=segmentTS(Y, k0=5)
#' # The transformed series z_t 
#' Z=Trans$X 
#' # Plot the cross correlogram of x_t and y_t
#' Z=data.frame(Z)
#' names(Z)=c("Z1","Z2","Z3","Z4","Z5","Z6")
#' # The cross correlogram of z_t shows 3-2-1 block pattern  
#' acfZ=acf(Z, plot=FALSE)
#' plot(acfZ, max.mfrow=6, xlab='', ylab='',  mar=c(1.8,1.3,1.6,0.5),
#'      oma=c(1,1.2,1.2,1), mgp=c(0.8,0.4,0),cex.main=1)
#' # Identify the permutation mechanism
#' permutation=permutationMax(Z) 
#' permutation$Groups 
#'
#'     
#' ## Example 2 (Example 6 of Chang et al.(2014)).
#' ## p=20, x_t consists of 5 independent subseries with 6, 5, 4, 3 and 2 components.    
#'
#' p=20;n=3000
#' # Generate x_t
#' X=mat.or.vec(p,n)
#' x=arima.sim(model=list(ar=c(0.5, 0.3), ma=c(-0.9, 0.3, 1.2,1.3)),n.start=500,n=n+5,sd=1)
#' for(i in 1:6) X[i,]=x[i:(n+i-1)]
#' x=arima.sim(model=list(ar=c(-0.4,0.5),ma=c(1,0.8,1.5,1.8)),n.start=500,n=n+4,sd=1)
#' for(i in 7:11) X[i,]=x[(i-6):(n+i-7)]
#' x=arima.sim(model=list(ar=c(0.85,-0.3),ma=c(1,0.5,1.2)), n.start=500,n=n+3,sd=1)
#' for(i in 12:15) X[i,]=x[(i-11):(n+i-12)]
#' x=arima.sim(model=list(ar=c(0.8,-0.5),ma=c(1,0.8,1.8)),n.start=500,n=n+2,sd=1)
#' for(i in 16:18) X[i,]=x[(i-15):(n+i-16)]
#' x=arima.sim(model=list(ar=c(-0.7, -0.5), ma=c(-1, -0.8)),n.start=500,n=n+1,sd=1)
#' for(i in 19:20) X[i,]=x[(i-18):(n+i-19)]
#' # Generate y_t 
#' A=matrix(runif(p*p, -3, 3), ncol=p)
#' Y=A%*%X  
#' Y=t(Y)
#' Trans=segmentTS(Y, k0=5)
#' # The transformed series z_t 
#' Z=Trans$X 
#' # Identify the permutation mechanism
#' permutation=permutationMax(Z) 
#' permutation$Groups  

 
## Output1: $NoGroups -- No. of groups with at least two components series
## Output2: $Nos_of_Members -- number of members in each of groups
## 				     with at least two members
## Output3: $Groups -- indices of components in each of groups with at least two members
## Output4: $maxcorr -- maximum correlation (over lags) of \eqn{p(p-1)/2} pairs in descending order
## Output5: $corrRatio -- ratios of succesive values from $maxCorr
## Output6: $NoConnectedPairs -- No. connected Pairs
# 
#
permutationMax <- function(X, Vol=FALSE, m=NULL) {
#
# X: nxp data matrix
# m: maximum lag used in calculating cross correlation coefficients
# 


## Step 0: prewhiten each columns of X
p=ncol(X)
n=nrow(X)
if(is.null(m))
{
m=10*log10(n/p)
}
if(!Vol)
{
R=5  
arOrder=rep(0, p)
 for(j in 1:p) { t=ar(X[,j],order.max=R)
                X[,j]=t$resid; arOrder[j]=t$order }
j=max(arOrder)
X=X[(j+1):n,]
}
if(Vol)
{
  ## Step 0: prewhiten each columns of X
  nanum=rep(0,p)
  for(j in 1:p) {options( warn = -1 ) 
                 t=garch(X[,j], order = c(1,1),trace=FALSE)
                  X[,j]=t$residuals
                  a=X[,j]
                  nanum[j]=length(a[is.na(X[,j])])    
                  
  }
  
  X=X[(max(nanum)+1):n,]
}
## Step 1: calculate max_k |\rho(k)| for each pair components

rho=acf(X,lag.max=m, plot=F) # rho$acf is an (m+1)xpxp array

p0=p*(p-1)/2 # total number of pairs of component series
M=vector(mode="numeric",length=p0) # max correlations between i-th and j-th component 
                         # over lags between -m to m, for 1 <= j < i <= p
for(i in 2:p) { for(j in 1:(i-1)) 
	M[(i-2)*(i-1)/2+j]=max(abs(rho$acf[,i,j]), abs(rho$acf[,j,i]))
	}
    # For a pxp matrix,  stack rows below the main diagoal together,
    # the (i,j)-th element, for i>j, is in the position (i-2)*(i-1)/2+j
# cat("STEP1","\n")

## Step 2: sorting p0 maximum correlation in descending order,
## find the ratio estimator for r
Ms=sort.int(M, decreasing=T, index.return=T)
    # Ms$x are sorted correlation, Ms$ix are the corresponding indices in M 
p1=as.integer(p0*0.75)
ratio=Ms$x[1:p1]/Ms$x[2:(p1+1)]
max0=max(ratio)
n=1:p1
r=n[ratio[n]==max0]; j=length(r); r=r[j]
# cat("STEP2","\n")

## Step 3: find the pairs corresponding to the r maximum max_k |\rho(k)|
h=mat.or.vec(p,1)
for(i in 2:p) h[i]=(i-2)*(i-1)/2
Inx=mat.or.vec(p,p); I=2:p
for(k in 1:r) { q=I[(Ms$ix[k]-h[I])>0]; s=length(q); i=q[s]; j=Ms$ix[k]-h[i]; Inx[i,j]=1}
# Now the entrices of Inx equal 1 are the positions with (i,j) connected,  
# and all other entrices are 0
# cat("STEP3","\n")

## Step 4: picking up the grouping from each columns of Inx, mark column with Index=1
##         with a group with at least two members, and Index=0 otherwise
G=mat.or.vec(p,p-1); 
Index=rep(0,p-1) 
N=rep(0,p-1)
# G[,j] records the components (from j-th column of Inx) to be grouped together with j
for(j in 1:(p-1)) { k=1
	for(i in (j+1):p) if(Inx[i,j]>0) { k=k+1; G[k,j]=i}
	if(k>1) { G[1,j]=j; Index[j]=1; N[j]=k }
}
# cat("STEP4","\n")

## Step 5: combining together any two groups obtained in Step 4 sharing
##         the same component
check=1
while(check>0) {  check=0
for(i in 1:(p-2)) { if(Index[i]==0) next
        for(j in (i+1):(p-1)) { if(Index[j]==0) next
                a=G[,i][G[,i]>0]; b=G[,j][G[,j]>0]
                c=c(a, b); d=length(c[duplicated(c)])
                if(d>0) {  # there are duplicated elements in a & b
                           a=unique(c) # picking up different elements from G[,i] & G[,j]
                           G[,i]=0; G[1:length(a),i]=sort(a); Index[j]=0
			   N[i]=length(a); N[j]=0;
			   check=1;
		}
        }
}
# cat("d=", d, "\n")
}
# cat("STEP5","\n")


## Step 6: Output
K=length(Index[Index==1])
if(K==0) stop("All component series are linearly independent")
Group=mat.or.vec(p,K+1)
k=1
for(j in 1:(p-1)) { if(Index[j]==1) { Group[,k]=G[,j]; k=k+1} }
cat("\n"); cat("No of groups with more than one members:", K, "\n")
cat("Nos of members in those groups:", N[N>0], "\n")
cat("No of connected pairs:", r, "\n")
#cat("Prewhited data are saved in the file Xpre.dat","\n\n")
#write.table(X, "Xpre.dat", row.names=F, col.names=F)
output=list(NoGroups=K, Nos_of_Members=N[N>0], Groups=Group[,1:K], maxcorr=Ms$x, corrRatio=ratio,NoConnectedPairs=r,Xpre=X)
return(output)
}
		      	
