#' Segment Multivariate Volatility Processes 
#'
#' Calculate linear transformation of the \eqn{p}-variate volatility processes {y_t} such that the transformed volatility process {x_t}=B{y_t} can be segmented into \eqn{q} lower-dimensional processes, and there exist no \emph{conditional} cross correlations across those \eqn{q} processes.
#'
##  When \eqn{p} is small, thresholding is not required. However, when \eqn{p} is large, it is necessary to use the thresholding method, see more information in Chang et al. (2014). 
#'
#' @param Y   a data matrix with \eqn{n} rows and \eqn{p} columns, where \eqn{n} is the sample size and \eqn{p} is the dimension of the time series.
#' @param k0  a positive integer specified to calculate Wy. 
 
#' @return An object of class "segmentVOL" is a list containing the following components:
#' \item{B}{the \eqn{p} by \eqn{p} transformation matrix such that {x_t}=B{y_t}}
#' \item{X}{the transformed series with \eqn{n} rows and \eqn{p} columns} 
## @note This is the first step to transform the time series. The second step is grouping the transformed time series, see \code{\link{permutationMax}}, \code{\link{permutationFDR}}.       
#' 
#' @author Jinyuan Chang, Bin Guo and Qiwei Yao
#' 
#' @references Chang, J., Guo, B. and Yao, Q. (2014).  \emph{Segmenting Multiple Time Series by Contemporaneous Linear Transformation: PCA for Time Series}. Available at \url{http://arxiv.org/abs/1410.2323}.
#' 
##             Cai, T. and Liu, W. (2011). \emph{Adaptive thresholding for sparse covariance matrix estimation}.  Journal of the American Statistical Association 106: 672-684.  
## @family aggregate functions
#' @seealso \code{\link{segmentTS}}      
#' @examples
#' ## Example 7 of Chang et al.(2014)
#' ## Segmenting the returns of the 6 stocks     
#' 
#' require(tseries)
#' data(returns)
#' Y=returns
#' n=dim(Y)[1]; p=dim(Y)[2]
## # The ACF plot of the residuals after prewhitening the original data by GARCH(1,1)
## nanum=rep(0,p)
## for(j in 1:p) { t=garch(Y_0[,j], order = c(1,1), trace=FALSE)
##               Y_0[,j]=t$residuals
##                a=Y_0[,j]
##               nanum[j]=length(a[is.na(Y_0[,j])]) }
## Y=Y_0[(max(nanum)+1):n,]
## colnames(Y)=c("Y1","Y2","Y3","Y4","Y5","Y6")
## t=acf(Y)
## plot(t, max.mfrow=6, xlab='', ylab='', mar=c(1.8,1.3,1.6,0.5),
##     oma=c(1,1.2,1.2,1), mgp=c(0.8,0.4,0), cex.main=1.0)
#' # Carry out the transformation procedure     
#' Trans=segmentVOL(Y,5)
## B=Trans$B    
#' X_0=data.frame(Trans$X)
#' X_1=X_0
#' # The ACF plot of the residuals after prewhitening the transformed data by GARCH(1,1)
#' nanum=rep(0,p)
#' for(j in 1:p) {options( warn = -1 ) 
#'                t=garch(X_1[,j], order = c(1,1),trace=FALSE)
#'                X_1[,j]=t$residuals
#'                a=X_1[,j]
#'                nanum[j]=length(a[is.na(X_1[,j])]) }
#' X=X_1[(max(nanum)+1):n,]
#' colnames(X)=c("X1","X2","X3","X4","X5","X6")
#' t=acf(X,plot=FALSE)
#' plot(t, max.mfrow=6, xlab='', ylab='',  mar=c(1.8,1.3,1.6,0.5),
#'     oma=c(1,1.2,1.2,1), mgp=c(0.8,0.4,0),cex.main=1.0,ylim=c(0,1))
#' # Identify the permutation mechanism
#' permutation=permutationMax(X_0,Vol=TRUE) 
#' permutation$Groups  
#' options( warn = 0) 

## This code is used to segment the multivariate volatility process
## Calculate segmentation transform X=BY via (i) stardardize y_t first
## and then the transformation in Step 1 of Chang, Guo and Yao (2014)

segmentVOL <- function(Y,k0) {

  # Part I -- standardize Y  such that var(y_t)=I_p
  
  n=nrow(Y)
  p=ncol(Y)
  
  Yacf=acf(Y, plot=F)
  M=var(Y) # M=var(y_t)
  t=eigen(M, symmetric=T)
  ev=sqrt(t$values) # square root of eigenvalues
  G=as.matrix(t$vectors)
  D=G*0; for(i in 1:p) { if(ev[i]>0) D[i,i]=1/ev[i]
                         else { print("Data are degenerate"); quit("yes")}
  }
  M1=G%*%D%*%t(G) # M1=var(y_t)^{-1/2}
  Y1=M1%*%t(Y)  
  Y=t(Y1) # Y is standardized now: var(y_t)=I_p
  
  Y=t(Y)

  # Part II -- Apply the transformation to recover x_t
  mean_y<-rowMeans(Y)
   
  
  #dyn.load("volsegment.dll")   # The following part is used to obtain W_y
  storage.mode(Y)<-"double"
  storage.mode(mean_y)<-"double"
  storage.mode(p)<-"integer"
  storage.mode(n)<-"integer"
  storage.mode(k0)<-"integer"
  res<-.Fortran("volsegment",Y,mean_y,n,p,k0,res=numeric(p^2))$res 

  W_y=matrix(res,nrow=p)

  t=eigen(W_y, symmetric=T)
  G=as.matrix(t$vectors)
  Y1=t(G)%*%Y 
  X=t(Y1) # segmented series
  B=t(G)%*%M1 # transformation matrix x_t = B y_t, does not include permutation in Step 2
  output=list(B=B, X=X)
}
