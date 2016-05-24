#' Segment Multivariate Time Series 
#'
#' Calculate linear transformation of the \eqn{p}-variate time series {y_t} such that the transformed series {x_t}=B{y_t} is segmented into several
#' lower-dimensional subseries, and those subseries are uncorrelated with
#' each other both contemporaneously and serially.
#'
#' When \eqn{p} is small, thresholding is not required. However, when \eqn{p} is large, it is necessary to use the thresholding method, see more information in Chang et al. (2014). 
#'
#' @param Y   a data matrix with \eqn{n} rows and \eqn{p} columns, where \eqn{n} is the sample size and \eqn{p} is the dimension of the time series.
#' @param k0  a positive integer specified to calculate Wy. See (2.5) in Chang et al. (2014).
#' @param thresh   logical. If \code{FALSE} (the default), no thresholding will be applied. If \code{TRUE}, a 
#'                 thresholding method will be applied first to estimate Wy, see (3.4) and (3.5) in Chang et al. (2014).     
#' @param tuning.vec  the value of thresholding parameter \eqn{\lambda}. The thresholding level is specified by             
#'                  \deqn{ u = \lambda {(log p/n)^(1/2)}.}
#'                    Default value is 2. If \code{tuning.vec} is a vector, then a cross validation method proposed in Cai and Liu (2011) will be used 
#'                    to choose the best tuning parameter. 
#'                     
#' @param K   the number of folders used in the cross validation, the default is 5. It is required when \code{thresh} is \code{TRUE}.                          
#' @return An object of class "segmentTS" is a list containing the following components:
#' \item{B}{the \eqn{p} by \eqn{p} transformation matrix such that {x_t}=B{y_t}}
#' \item{X}{the transformed series with \eqn{n} rows and \eqn{p} columns} 
#' @note This is the first step to transform the time series. The second step is grouping the transformed time series, see \code{\link{permutationMax}}, \code{\link{permutationFDR}}.       
#' 
#' @author Jinyuan Chang, Bin Guo and Qiwei Yao
#' 
#' @references Chang, J., Guo, B. and Yao, Q. (2014).  \emph{Segmenting Multiple Time Series by Contemporaneous Linear Transformation: PCA for Time Series}. Available at \url{http://arxiv.org/abs/1410.2323}.
#' 
#'             Cai, T. and Liu, W. (2011). \emph{Adaptive thresholding for sparse covariance matrix estimation}.  Journal of the American Statistical Association 106: 672-684.  
## @family aggregate functions
#' @seealso \code{\link{permutationMax}}, \code{\link{permutationFDR}}.       
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
#' # Plot the cross correlogram of z_t and y_t
#' Y=data.frame(Y);Z=data.frame(Z)
#' names(Y)=c("Y1","Y2","Y3","Y4","Y5","Y6")
#' names(Z)=c("Z1","Z2","Z3","Z4","Z5","Z6")
#' # The cross correlogram of y_t shows no block pattern 
#' acfY=acf(Y) 
#' # The cross correlogram of z_t shows 3-2-1 block pattern  
#' acfZ=acf(Z)
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
#' # Plot the cross correlogram of x_t and y_t
#' Y=data.frame(Y);Z=data.frame(Z)
#' namesY=NULL;namesZ=NULL
#' for(i in 1:p)
#' {
#'    namesY=c(namesY,paste0("Y",i))
#'    namesZ=c(namesZ,paste0("Z",i))
#' }  
#' names(Y)=namesY;names(Z)=namesZ
#' # The cross correlogram of y_t shows no block pattern 
#' acfY=acf(Y, plot=FALSE)
#' plot(acfY, max.mfrow=6, xlab='', ylab='',  mar=c(1.8,1.3,1.6,0.5), 
#'      oma=c(1,1.2,1.2,1), mgp=c(0.8,0.4,0),cex.main=1)    
#' # The cross correlogram of z_t shows 6-5-4-3-2 block pattern  
#' acfZ=acf(Z, plot=FALSE)
#' plot(acfZ, max.mfrow=6, xlab='', ylab='',  mar=c(1.8,1.3,1.6,0.5),
#'      oma=c(1,1.2,1.2,1), mgp=c(0.8,0.4,0),cex.main=1)
#' # Identify the permutation mechanism
#' permutation=permutationMax(Z) 
#' permutation$Groups  




segmentTS <- function(Y,k0,thresh=FALSE,tuning.vec=2,K=5)
{
  # Part I -- standardize Y  such that var(y_t)=I_p
  n=nrow(Y)
  p=ncol(Y)
  # M=var(y_t)
  M=var(Y)            
  t=eigen(M, symmetric=T)
  # square root of eigenvalues
  ev=sqrt(t$values)   
  G=as.matrix(t$vectors)
  D=G*0; for(i in 1:p) 
  { 
    if(ev[i]>0) D[i,i]=1/ev[i]
    else { print("Data are degenerate"); quit("yes")}
  }
  Y0=Y
  # M1=var(y_t)^{-1/2}
  M1=G%*%D%*%t(G)    
  Y1=M1%*%t(Y)  
  # Y is standardized now: var(y_t)=I_p
  Y=Y1                
  
  # Part II -- Apply the transformation to recover x_t
  
  Wy=diag(rep(1,p))
  mean_y<-rowMeans(Y)
  storage.mode(Y)<-"double"
  storage.mode(mean_y)<-"double"
  storage.mode(p)<-"integer"
  storage.mode(n)<-"integer"
  
  # Segment without thresholding  
  
  if(!thresh)          
  {

    for(k in 1:k0) {
      S=cov(t(Y[,1:(n-k)]),t(Y[,(1+k):n])); Wy=Wy+S%*%t(S)
    } 
  }
  
  # Segment with thresholding 
  
  if(thresh)            
  { 
     if(length(tuning.vec)<1)
     {
       print("Need to specify the thresholding parameter"); quit("yes")
     }
     if(length(tuning.vec)==1)
     {
       K=1               # Using the fixed thresholding 
     }
     
     # Using the cross validation for thresholding 
     
     for(k in 1:k0)         
     { 
       error=NULL
       storage.mode(k)<-"integer"
       
       # To select proper threshold parameter
       
       for(v in 1:K)
       {
         sample1=sample(1:n,size=n/2)
         sample2=c(1:n)[-sample1]
         sampleY1=Y[,sample1]
         sampleY2=Y[,sample2]
         
         mean_y1<-rowMeans(sampleY1)
         mean_y2<-rowMeans(sampleY2)
         
         storage.mode(mean_y1)<-"double"
         storage.mode(sampleY1)<-"double"
         storage.mode(mean_y2)<-"double"
         storage.mode(sampleY2)<-"double"
         n1=ceiling(n/2)
         storage.mode(n1)<-"integer"
         storage.mode(p)<-"integer"
         
         errors=NULL
         
         for(d in 1:length(tuning.vec))
         { 
           delta1=tuning.vec[d]
           storage.mode(delta1)<-"double"
           
           res1<-.Fortran("segment",sampleY1,mean_y1,k,n1,p,res=numeric(p^2))$res 
           Sigma_y1<-matrix(res1,nrow=p,byrow =TRUE)
           
           res2<-.Fortran("segment",sampleY2,mean_y2,k,n1,p,res=numeric(p^2))$res 
           Sigma_y2<-matrix(res2,nrow=p,byrow =TRUE)
           
           storage.mode(Sigma_y1)<-"double"
           res<-.Fortran("thresh",Sigma_y1,sampleY1,mean_y1,k,n1,p,delta1,res1=numeric(p^2))$res1
           Sigma_ythres1<-matrix(res,nrow=p,byrow =TRUE)
           
           errors=c(errors,(norm((Sigma_ythres1-Sigma_y2),type="F"))^2)   
           }  
         error=rbind(error,errors)
       }
       errormean=colMeans(error)
       d=which.min(errormean) 
       deltafinal=tuning.vec[d]         
       
       # Find the best tuning parameter 
       
       res<-.Fortran("segment",Y,mean_y,k,n,p,res=numeric(p^2))$res 
       Sigma_y<-matrix(res,nrow=p,byrow =TRUE)
       storage.mode(Sigma_y)<-"double"
       storage.mode(deltafinal)<-"double"
       
       # Carry out the final thresholding 
       
       res<-.Fortran("thresh",Sigma_y,Y,mean_y,k,n,p,deltafinal,res1=numeric(p^2))$res1
       Sigma_ynew<-matrix(res,nrow=p,byrow =TRUE)
       Wy=Wy+Sigma_ynew%*%t(Sigma_ynew)  
     }
     
     
  }
  t=eigen(Wy, symmetric=T)
  G=as.matrix(t$vectors)
  Y1=t(G)%*%Y 
  # segmented series
  X=t(Y1)  
  # transformation matrix x_t = B y_t, does not include permutation in Step 2
  B=t(G)%*%M1       
   
  output=list(B=B, X=X)
  return(output)
  
}
