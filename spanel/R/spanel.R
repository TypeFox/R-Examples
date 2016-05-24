#' The contiguity matrix
#'
#' @docType data
#' @name  usaww
#' @usage data(usaww)
NULL
#' US States Production
#'
#' \itemize{
#' \item{state}{the state}
#' \item{year}{the year}
#' \item{pcap}{private capital stock}
#' \item{hwy}{highway and streets}
#' \item{water}{water and sewer facilities}
#' \item{util}{other public buildings and structures}
#' \item{pc}{public capital}
#' \item{gsp}{gross state products}
#' \item{emp}{labor input measured by the employement in non--agricultural payrolls}
#' \item{unemp}{state unemployment rate}
#' }
#'
#' @docType data
#' @name Produc
#' @usage data(Produc)
#' @format A data frame with 816 rows and 10 variables
NULL


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
spfe<-function(y,x,w,n,t){
  x<-x
  y<-y
  ww<-as.matrix(w)
  h<-kronecker(ww,diag(t))
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
  H<-cbind(x,h%*%x,h^2%*%x)
  Z<-cbind(x,h%*%y)
  zfe<-Q%*%Z
  hfe<-Q%*%H
  yfe<-Q%*%y
  th<-crossprod(hfe)
  thz<-crossprod(hfe,zfe)
  zfehat<-hfe%*%ginv(th)%*%thz
  beta<-ginv(crossprod(zfehat,zfe))%*%crossprod(zfehat,yfe)
  ieffct<-P%*%(y-Z%*%beta)
  yhat<-Z%*%beta+ieffct
  #nn<-ncol(Z)-2
  #s<-length(beta)-2
  #j<-cbind(Z[,1:nn],Z[,ncol(Z)])
  #jj<-rbind(beta[1:s],beta[length(beta)])
  #yhat<-ginv(diag(N)-beta[length(beta)-1]*h)%*%(jj%*%j+ieffect)
  #----------------
  res<-y-yhat
  resvar<-sum(res^2)/(N-k-n+1)
  varcoef<-resvar*ginv(crossprod(zfehat,zfe))
  scr<-sum(res^2)
  sct<-crossprod(yfe,yfe)
  rsq<-(cor(yhat,yfe))^2
  std<-sqrt(diag(varcoef[-1,-1]))
  tval<-beta[-1]/std
  pval = 2*pt(-abs(tval), df=df)
  list(coefficients=beta[-1],resvar=resvar,vcov=varcoef,std=std,tval=tval,pval=pval,Rsq=rsq,df=df,n=n,t=t,N=N,k=k)
}


spbe<-function(y,x,w,n,t){
  x<-x
  y<-y
  w<-as.matrix(w)
  h<-kronecker(w,diag(t))
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
  H<-cbind(x,h%*%x,h^2%*%x)
  Z<-cbind(x,h%*%y)
  zbe<-P%*%Z
  ybe<-P%*%y
  hbe<-P%*%H
  zbehat<-hbe%*%ginv(crossprod(hbe))%*%crossprod(hbe,zbe)
  beta<-ginv(crossprod(zbehat,zbe))%*%crossprod(zbehat,ybe)
  yhat<-zbe%*%beta
  res<-y-yhat
  resvar<-sum(res^2)/(n-k)
  varcoef<-resvar*ginv(crossprod(zbehat,zbe))
  scr<-sum(res^2)
  sct<-crossprod(ybe,ybe)
  rsq<-(cor(yhat,ybe))^2
  std<-sqrt(diag(varcoef))
  tval<-beta/std
  pval = 2*pt(-abs(tval), df=df)
  list(coefficients=beta,resvar=resvar,vcov=varcoef,std=std,tval=tval,pval=pval,Rsq=rsq,df=df,n=n,t=t,N=N,k=k)
}

sprandom<-function(y,x,w,n,t){
  N<-n*t
  k<-ncol(x)
  w<-as.matrix(w)
  h<-kronecker(w,diag(t))
  df<-nrow(x)-ncol(x)
  ones<-matrix(1,t,t)
  jbar<-ones/t
  E<-diag(t)-jbar
  P<-kronecker(diag(n),jbar)
  fx<-spfe(y,x,w,n,t)
  be<-spbe(y,x,w,n,t)
  sigmaf<-fx$resvar
  sigmab<-be$resvar
  sigmamu<-(sigmab-sigmaf)/t
  rhomu<-sigmamu/(sigmamu+sigmaf)
  theta<-1-sqrt(sigmaf/(t*sigmamu+sigmaf))
  Qre<-diag(N)-theta*P
  H<-cbind(x,h%*%x,h^2%*%x)
  Z<-cbind(x,h%*%y)
  yre<-Qre%*%y
  xre<-Qre%*%Z
  hre<-Qre%*%H
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
#' data(Produc)
#' data("usaww")
#' #fit the fixed function
#' fx<-span(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp,Produc,usaww,n=48,t=17,model="fe")
#' # fit the random function
#'  ran<-span(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp,Produc,usaww,n=48,t=17,model="re")
#' # the Hausman test
#' hausman(fx,ran)
#' @export
hausman<-function(fixed,random){
  fecoef<-fixed$coefficients
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
span<- function(x,...){UseMethod("span") }



span.default <- function(y,x,w,n,t,model=c('fe','be','re'),...)
{

  t<-t
  n<-n
  x<-as.matrix(x)
  h<-as.matrix(w)
  y<-as.numeric(y)
  if(model=="fe"){
    est <- spfe(y,x,w,n,t)
  }
  if(model=="be"){
    est <- spbe(y,x,w,n,t)
  }
  if(model=="re"){
    est <- sprandom(y,x,h,n,t)
  }
  est$call <- match.call()
  class(est) <- "span"
  est


}
print.span <- function(x,...)
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
summary.span<-function(object,...)
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
#' @param formula log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp
#' @param w is the contiguity matrix
#' @param data the dataframe
#' @param n the number of section
#' @param t the time per section
#' @param model  "fe" for fixed effect "be" for between and "re" for random effect
#' @param ... not used
#' @export
span.formula <-function(formula,data=list(),w,n,t,model=c("fe","be","re"),...)
{
  t<-t
  n<-n
  mf <- model.frame(formula=formula, data=data)
  x <- model.matrix(attr(mf, "terms"), data=mf)
  y <- model.response(mf)
  h<-w
  if(model=="fe"){
    est <- span.default(y,x,h,n,t,model="fe",...)
  }
  if(model=="be"){
    est <- span.default(y,x,h,n,t,model="be",...)
  }
  if(model=="re"){
    est <- span.default(y,x,h,n,t,model="re",...)
  }
  est$call <- match.call()
  est$equa <- formula
  est
}
