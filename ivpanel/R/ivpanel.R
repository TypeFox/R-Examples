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
ivfixed<-function(y,x,h,n,t){
  x<-x
  y<-y
  h<-h
  t<-t
  n<-n
  N<-n*t
  k<-ncol(x)
  df<-nrow(x)-ncol(x)
  ones<-matrix(1,t,t)
  jbar<-ones/t
  E<-diag(t)-jbar
  P<-kronecker(diag(n),jbar)
  Q<-diag(N)-P
  xfe<-Q%*%x
  yfe<-Q%*%y
  hfe<-Q%*%h
  txfe<-crossprod(xfe,xfe)
  tyfe<-crossprod(xfe,yfe)
  thfe<-crossprod(hfe,hfe)
  txhfe<-crossprod(hfe,xfe)
  xfehat<-hfe%*%ginv(thfe)%*%txhfe
  beta<-ginv(crossprod(xfehat,xfe))%*%crossprod(xfehat,yfe)
  ieffect<-P%*%(y-x%*%beta)
  yhat<-x%*%beta+ieffect
  res<-y-yhat
  resvar<-sum(res^2)/(n*t-k-n+1)
  varcoef<-resvar*ginv(crossprod(xfehat,xfe))
  scr<-sum(res^2)
  sct<-crossprod(yfe,yfe)
  rsq<-(cor(yhat,yfe))^2
  std<-sqrt(diag(varcoef))
  tval<-beta/std
  pval = 2*pt(-abs(tval), df=df)
  list(coefficients=beta,resvar=resvar,vcov=varcoef,std=std,tval=tval,pval=pval,Rsq=rsq,df=df,n=n,t=t,N=N,k=k)
}


ivbe<-function(y,x,h,n,t){
  x<-x
  y<-y
  h<-h
  t<-t
  n<-n
  N<-n*t
  k<-ncol(x)
  df<-nrow(x)-ncol(x)
  ones<-matrix(1,t,t)
  jbar<-ones/t
  E<-diag(t)-jbar
  P<-kronecker(diag(n),jbar)
  Q<-diag(N)-P
  xfe<-P%*%x
  yfe<-P%*%y
  hfe<-P%*%h
  txfe<-crossprod(xfe,xfe)
  tyfe<-crossprod(xfe,yfe)
  thfe<-crossprod(hfe,hfe)
  txhfe<-crossprod(hfe,xfe)
  xfehat<-hfe%*%ginv(thfe)%*%txhfe
  beta<-ginv(crossprod(xfehat,xfe))%*%crossprod(xfehat,yfe)
  yhat<-xfe%*%beta
  res<-y-yhat
  resvar<-sum(res^2)/(n-k)
  varcoef<-resvar*ginv(crossprod(xfehat,xfe))
  scr<-sum(res^2)
  sct<-crossprod(yfe,yfe)
  rsq<-(cor(yhat,yfe))^2
  std<-sqrt(diag(varcoef))
  tval<-beta/std
  pval = 2*pt(-abs(tval), df=df)
  list(coefficients=beta,resvar=resvar,vcov=varcoef,std=std,tval=tval,pval=pval,Rsq=rsq,df=df,n=n,t=t,N=N,k=k)
}

ivrandom<-function(y,x,h,n,t){
  N<-n*t
  k<-ncol(x)
  df<-nrow(x)-ncol(x)
  ones<-matrix(1,t,t)
  jbar<-ones/t
  E<-diag(t)-jbar
  P<-kronecker(diag(n),jbar)
 fx<-ivfixed(y,x,h,n,t)
 be<-ivbe(y,x,h,n,t)
  sigmaf<-fx$resvar
  sigmab<-be$resvar
  sigmamu<-(sigmab-sigmaf)/t
  rhomu<-sigmamu/(sigmamu+sigmaf)
  theta<-1-sqrt(sigmaf/(t*sigmamu+sigmaf))
  Qre<-diag(N)-theta*P
  yre<-Qre%*%y
  xre<-Qre%*%x
  hre<-Qre%*%h
  kk<-ncol(xre)
  xrehat<-hre%*%(ginv(crossprod(hre))%*%crossprod(hre,xre))
  rbeta<- ginv(crossprod(xrehat,xre))%*%crossprod(xrehat,yre)
  ryhat<-xre%*%rbeta
  rres<-yre-ryhat
  rresvar<-sum(rres^2)/(N-kk)
  rvarcoef<-rresvar*ginv(crossprod(xrehat,xre))
  rrsq<-(cor(ryhat,yre))^2
  rstd<-sqrt(diag(rvarcoef))
  rtval<-rbeta/rstd
  rpval = 2*pt(-abs(rtval), df=df)
  list(coefficients=rbeta,vcov=rvarcoef,std=rstd,tval=rtval,pval=rpval,Rsq=rrsq,df=df,n=n,t=t,N=N,k=k)
}

#' Hausman test
#' 
#' @param fixed is the fixed effect object function
#' @param random is the random effect object function
#' @return Chisq the hausman statistic
#' @return P-value the probability value
#' @return df the degree of freedom
#' @examples
#' pib<-as.matrix(c(12,3,4,0.4,0.7,5,0.7,0.3,0.6,89,7,8,45,7,4,5,0.5,5),nrows=18,ncols=1)
#' tir<-as.matrix(c(12,0.3,4,0.4,7,12,3.0,6.0,45,7.0,0.8,44,65,23,4,6,76,9),nrows=18,ncols=1)
#' inf<-as.matrix(c(1.2,3.6,44,1.4,0.78,54,0.34,0.66,12,0.7,8.0,12,65,43,5,76,65,8),nrows=18,ncols=1)
#' npl<-as.matrix(c(0.2,3.8,14,2.4,1.7,43,0.2,0.5,23,7.8,88,36,65,3,44,65,7,34),nrows=18,ncols=1)
#' #create a data frame 
#' mdata<-data.frame(p=pib,t=tir,int=inf,np=npl)
#' #fit the fixed function  
#' fx<-ivpan(t~p+int|p+np,mdata,n=6,t=3,model="fe")
#' # fit the random function
#' ran<-ivpan(t~p+int|p+np,mdata,n=6,t=3,model="re")
#' # the Hausman test
#' hausman(fx,ran)
#' @export
hausman<-function(fixed,random){
  fecoef<-fixed$coefficients[-1]
  recoef<-random$coefficients[-1]
  fcoefvar<-fixed$vcov[-1,-1]
  rcoefvar<-random$vcov[-1,-1]
  h<-t(fecoef-recoef)%*%ginv(fcoefvar-rcoefvar)%*%(fecoef-recoef)
  df=length(fecoef)
  pval<- 1-pchisq(h,df=df)
  cat("\nHausman Test")
  cat("\nChisq:",h,"\nP-value:",pval,"\ndf:",df)
}


#' method
#' 
#' @author Zaghdoudi Taha
#' @param x a numeric design matrix for the model.
#' @param ... not used
#' @export
ivpan<- function(x,...){UseMethod("ivpan") }



ivpan.default <- function(y,x,h,n,t,model=c('fe','be','re'),...)
{
  
  t<-t
  n<-n
  x<-as.matrix(x)
  h<-as.matrix(h)
  y<-as.numeric(y)
  if(model=="fe"){
    est <- ivfixed(y,x,h,n,t)
  }
  if(model=="be"){
    est <- ivbe(y,x,h,n,t)
  }
  if(model=="re"){
    est <- ivrandom(y,x,h,n,t)
  }
  est$call <- match.call()
  class(est) <- "ivpan"
  est
  
  
}
print.ivpan <- function(x,...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
}

#' Summary
#' 
#' @param object is the object of the function
#' @param ... not used
#' @export
summary.ivpan<-function(object,...)
{
  res <- cbind(object$coefficients,object$std, object$tval,object$pval )
  colnames(res) <- c("Estimates", "Std.Err", "T-value", "P-Value")
  cat("Formula:")
  print(object$equa)
  cat("\nBalanced Panel:","n:",object$n,"t:",object$t,"N:",object$N,"\n")
  cat("Rsquared :",object$Rsq,"\n")
  printCoefmat(res,has.Pvalue=TRUE)
 
    
}
#' formula
#' 
#' @param formula PIB~INF+TIR|Cap+m2r "|" rhs is the instrumental variables
#' @param data the dataframe 
#' @param n the number of section
#' @param t the time per section
#' @param model  "fe" for fixed effect "be" for between and "re" for random effect
#' @param ... not used
#' @import Formula
#' @export
ivpan.formula <-function(formula,data=list(),n,t,model=c("fe","be","re"),...)
{
  t<-t
  n<-n
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0)
  mf <- mf[c(1, m)]
  f <- Formula(formula)
  mf[[1]] <- as.name("model.frame")
  mf$formula <- f
  mf <- eval(mf, parent.frame())
  x <- model.matrix(f, data = mf, rhs = 1)
  h <- model.matrix(f, data = mf, rhs = 2)
  y <- model.response(mf)
  if(model=="fe"){
    est <- ivpan.default(y,x,h,n,t,model="fe",...)
  }
  if(model=="be"){
    est <- ivpan.default(y,x,h,n,t,model="be",...)
  }
  if(model=="re"){
    est <- ivpan.default(y,x,h,n,t,model="re",...)
  }
  est$call <- match.call()
  est$equa <- formula
  est
}