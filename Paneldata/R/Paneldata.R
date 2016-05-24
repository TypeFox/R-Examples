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
#' fixed effect function
#' 
#' @author Zaghdoudi Taha
#' @param y vector of dependent variable
#' @param x matrix of independents variables
#' @param n number of sections
#' @param t times per section
#' @return Coefficients a named vector of coefficients
#' @return vcov covariance matrix of coefficients
#' @return df the degree of freedom
#' @export
#' @examples
#' pib<-as.matrix(c(12,3,4,0.4,0.7,5,0.7,0.3,0.6,89,7,8,45,7,4,5,0.5,5),nrows=18,ncols=1)
#' tir<-as.matrix(c(12,0.3,4,0.4,7,12,3.0,6.0,45,7.0,0.8,44,65,23,4,6,76,9),nrows=18,ncols=1)
#' inf<-as.matrix(c(1.2,3.6,44,1.4,0.78,54,0.34,0.66,12,0.7,8.0,12,65,43,5,76,65,8),nrows=18,ncols=1)
#' npl<-as.matrix(c(0.2,3.8,14,2.4,1.7,43,0.2,0.5,23,7.8,88,36,65,3,44,65,7,34),nrows=18,ncols=1)
#' # create data frame
#' mdata<-data.frame(p=pib,ti=tir,int=inf,np=npl)
#' # create the designed matrix for the model
#' d<-matrix(c(mdata$p,mdata$int,mdata$np),18, 3)
#' # Fit a fixed model
#' fx<-fixed(mdata$p,d,n=6,t=3)
#' fx
fixed<-function(y,x,n,t){
  x<-x
  y<-y
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
  txfe<-crossprod(xfe,xfe)
  tyfe<-crossprod(xfe,yfe)
  beta<-ginv(txfe)%*%tyfe
  ieffect<-P%*%(y-x%*%beta)
  yhat<-x%*%beta+ieffect
  res<-y-yhat
  resvar<-sum(res^2)/(n*(t-1)-k+1)
  varcoef<-resvar*ginv(txfe)
  list(coefficients=beta,vcov=varcoef,df=df)
}
#' random effect function
#' 
#' @author Zaghdoudi Taha
#' @param y vector of dependent variable
#' @param x matrix of independents variables
#' @param n number of sections
#' @param t times per section
#' @return Coefficients: a named vector of coefficients
#' @return vcov: covariance matrix of coefficients
#' @return std: a named vector of standard errors 
#' @return stats: a named vector of students statistics
#' @export
#' @examples
#' pib<-as.matrix(c(12,3,4,0.4,0.7,5,0.7,0.3,0.6,89,7,8,45,7,4,5,0.5,5),nrows=18,ncols=1)
#' tir<-as.matrix(c(12,0.3,4,0.4,7,12,3.0,6.0,45,7.0,0.8,44,65,23,4,6,76,9),nrows=18,ncols=1)
#' inf<-as.matrix(c(1.2,3.6,44,1.4,0.78,54,0.34,0.66,12,0.7,8.0,12,65,43,5,76,65,8),nrows=18,ncols=1)
#' npl<-as.matrix(c(0.2,3.8,14,2.4,1.7,43,0.2,0.5,23,7.8,88,36,65,3,44,65,7,34),nrows=18,ncols=1)
#' mdata<-data.frame(p=pib,t=tir,int=inf,np=npl)
#' # create the designed matrix for the model
#' d<-matrix(c(mdata$p,mdata$int,mdata$np),18, 3)
#' # Fit a random model
#' rx<-Rand(mdata$p,d,n=6,t=3)
#' rx
Rand<-function(y,x,n,t){
  x<-x
  y<-y
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
  txfe<-crossprod(xfe,xfe)
  tyfe<-crossprod(xfe,yfe)
  beta<-ginv(txfe)%*%tyfe
  ieffect<-P%*%(y-x%*%beta)
  yhat<-x%*%beta+ieffect
  res<-y-yhat
  resvar<-sum(res^2)/(n*(t-1)-k+1)
  #between effect
  xbe<-P%*%x
  ybe<-P%*%y
  txbe<-crossprod(xbe,xbe)
  tybe<-crossprod(xbe,ybe)
  betabe<-ginv(txbe)%*%tybe
  ieffectbe<-P%*%(y-x%*%betabe)
  yhatbe<-x%*%betabe+ieffectbe
  resbe<-ybe-yhatbe
  resvarbe<-sum(resbe^2)/(n-k)
  varcoefbe<-resvarbe*ginv(txbe)
  stderrbe<-sqrt(diag(varcoefbe))
  statbe<-betabe/stderrbe  
  # random effect
  sigmau<-resvarbe-(resvar/t)
  o<-resvar/(t*sigmau+resvar)
  phi<-sqrt(o)
  teta<-1-phi
  Qre<-diag(N)-teta*P
  yre<-Qre%*%y
  xre<-Qre%*%x
  txre<-crossprod(xre,xre)
  tyre<-crossprod(xre,yre)
  betare<-ginv(txre)%*%tyre
  yhatre<-xre%*%betare
  resre<-yre-yhatre
  resvarre<-sum(resre^2)/(N-k)
  varcoefre<-resvarre*ginv(txre)
  stderrre<-sqrt(diag(varcoefre))
  statre<-betare/stderrre
  list(coefficients=betare,Std=stderrre,stat=statre,vcov=varcoefre,df=df)
  
  
}
#' method
#' 
#' @author Zaghdoudi Taha
#' @param x a numeric design matrix for the model.
#' @param ... not used
#' @export
Paneldata <- function(x,...){UseMethod("Paneldata") }



Paneldata.default <- function(y,x,n,t,model=c('fe','re'),...)
{
  
  
  t<-t
  n<-n
  x<-as.matrix(x)
  y<-as.numeric(y)
  if(model=="fe"){
    est <- fixed(y,x,n,t)
  }
  if(model=="re"){
    est <- Rand(y,x,n,t)
  }
  est$fitted.values <- as.vector(x %*% est$coefficients)
  est$residuals <- y - est$fitted.values
  est$call <- match.call()
  class(est) <- "Paneldata"
  est
  
  
}
print.Paneldata <- function(x,...)
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
summary.Paneldata<-function(object,...)
{
  se <- sqrt(diag(object$vcov))
  tval <-coef(object)/se
  TAB <- cbind(Estimate = coef(object),StdErr = se, t.value = tval, p.value = 2*pt(-abs(tval), df=object$df))
  colnames(TAB) <- c("Estimate", "Std.Err", "T value", "Pr(>z)")
  res <- list(call=object$call,
              coefficients=TAB)
  
  res
}
#' formula
#' 
#' @param formula PIB~INF+TIR
#' @param data the dataframe 
#' @param n the number of section
#' @param t the time per section
#' @param model if fixed==>"fe" if randon==>"re"
#' @param ... not used
#' @export
Paneldata.formula <-function(formula, data=list(),n,t,model=c('fe','re'),...)
{
  t<-t
  n<-n
  mf <- model.frame(formula=formula, data=data)
  x <- model.matrix(attr(mf, "terms"), data=mf)
  y <- model.response(mf)
  if(model=="fe")
    est <- Paneldata.default(y,x,n,t,model="fe",...)
  if(model=="re")
    est <- Paneldata.default(y,x,n,t,model="re",...)
  est$call <- match.call()
  est$formula <- formula
  est
}