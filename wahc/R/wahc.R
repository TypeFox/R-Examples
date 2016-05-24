ginv <- function(x, tol = sqrt(.Machine$double.eps))
{
  ## Generalized Inverse of a Matrix
  dnx <- dimnames(x)
  if(is.null(dnx)) dnx <- vector("list", 2)
  s <- svd(x)
  nz <- s$d > tol * s$d[1]
  structure(
    if(any(nz)) s$v[, nz] %*% (t(s$u[, nz])/s$d[nz]) else x,
    dimnames = dnx[2:1])
}
pvmat<-function(x,n,t){
  k<-ncol(x)
  pv<-matrix(0,n*t,k)
  ones<-matrix(1,t,1)
  for(i in 1:n){
    s=(i-1)*t+1
    e=t*i
    pv[s:e,]<-kronecker(ones,t(colMeans(x[s:e,])))
  }
  return(pv)
}
ypvmat<-function(x,n,t){
  k<-ncol(x)
  pv<-matrix(0,n*t,1)
  ones<-matrix(1,t,1)
  for(i in 1:n){
    s=(i-1)*t+1
    e=t*i
    pv[s:e]<-kronecker(ones,t(mean(x[s:e])))
  }
  return(pv)
}
yqvmat<-function(x,n,t){
  y<-x-ypvmat(x,n,t)
  return(y)
}

qvmat<-function(x,n,t){
  y<-x-pvmat(x,n,t)
  return(y)
}
reshape<-function(x,n,t){
  xx<-matrix(x,n,t)
  y<-t(xx)
  return(y)
}
wahc.fit<-function(y,x,n,t){
 
  qvx<-qvmat(x,n,t)
  qvy<-yqvmat(y,n,t)
  b<-ginv(crossprod(qvx))%*%crossprod(qvx,qvy)
  e<-qvy-qvx%*%b
  bb<-matrix(0,n*t,ncol(qvx))
  for(i in 1: n){
    s=(i-1)*t+1
    ee=t*i
    bb[s:ee,]<-e[s:ee]%*%t(e[s:ee])%*%qvx[s:ee,]
  }
  sb<-crossprod(qvx,bb)
  varcoef<-ginv(crossprod(qvx))%*%sb%*%ginv(crossprod(qvx))
  sd<-sqrt(diag(varcoef))
  st<-b/sd
  df=nrow(qvx)-ncol(qvx)
  N=n*t
  pval<-2*pt(-abs(st), df=df)
  rsq<-(cor(qvx%*%b,qvy))^2
  res<-list(coefficients=b,std=sd,tstat=st,pv=pval,rsq=rsq,varcoef=varcoef,df=df,n=n,t=t,N=N)
  class(res)<- "whc"
  return(res)
}


#' Summary
#' 
#' @param object is the object of the function
#' @param ... not used
#' @export
summary.whc<-function(object,...)
{
  res<-cbind(object$coefficients,object$std,object$tstat,object$pv)
  colnames(res)<-c("Estimates","stderr","t-value","pvalue")
  cat("Formula:")
  print(object$formula)
  cat("\nBalanced Panel:","n:",object$n,"t:",object$t,"N:",object$N,"\nRsquared :", object$rsq,"\n")
  printCoefmat(res,has.Pvalue=TRUE)
}

#' Fitting the fixed effect panel data model with heteroskedasticity and autocorrelation correction
#' 
#' @usage whc(formula, data,n,t,...)
#' @param formula an object of class \code{\link{formula}} 
#' @param data the dataframe 
#' @param n the number of section
#' @param t the time per section
#' @param ... not used
#' @examples
#' # Create data
#' pib<-as.matrix(c(12,3,4,0.4,0.7,5,0.7,0.3,0.6,89,7,8,45,7,4,5,0.5,5),nrows=18,ncols=1)
#' tir<-as.matrix(c(12,0.3,4,0.4,7,12,3.0,6.0,45,7.0,0.8,44,65,23,4,6,76,9),nrows=18,ncols=1)
#' inf<-as.matrix(c(1.2,3.6,44,1.4,0.78,54,0.34,0.66,12,0.7,8.0,12,65,43,5,76,65,8),nrows=18,ncols=1)
#' npl<-as.matrix(c(0.2,3.8,14,2.4,1.7,43,0.2,0.5,23,7.8,88,36,65,3,44,65,7,34),nrows=18,ncols=1)
#' #create a data frame 
#' mdata<-data.frame(p=pib,t=tir,int=inf,np=npl)
#' #fit the model  
#' fx<-whc(p~int+t,mdata,n=6,t=3)
#' summary(fx)
#' @export
whc<-function(formula,data,n,t,...){
  n<-n
  t<-t
  mf <- model.frame(formula=formula, data=data)
  x <- model.matrix(attr(mf, "terms"), data=mf)
  y <- model.response(mf)
  est<- wahc.fit(y,x,n,t)
  est$formula<-formula
  est
}
