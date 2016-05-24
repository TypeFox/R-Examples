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
#' least square dummy variable
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
#' # Fit a least square dummy variable model
#' ls<-lsdv(mdata$p,d,n=6,t=3)
#' ls
lsdv<-function(y,x,n,t){
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
  scr<-sum(res^2)
  sct<-crossprod(yfe,yfe)
  rsq<-1-scr/sct
  #MCO normale
  xq<-qr(x)
  alfa<-solve.qr(xq,y)
  error<-y-x%*%alfa
  scre<-sum((error)^2)
  effet<-(scre-scr)*(n*t-n-k)/(scr*(n-1))
  dff<-n-1
  Dff<-n*t-n-k
  peffet<-pf(effet,dff,Dff, lower.tail = TRUE, log.p = FALSE)
  list(coefficients=beta,vcov=varcoef,Rsq=rsq,df=df,Effect=effet,Dff=Dff,Peffect=peffet)
}

  

#' method
#' 
#' @author Zaghdoudi Taha
#' @param x a numeric design matrix for the model.
#' @param ... not used
#' @export
Lsdv <- function(x,...){UseMethod("Lsdv") }



Lsdv.default <- function(y,x,n,t,...)
{
   
  t<-t
  n<-n
  x<-as.matrix(x)
  y<-as.numeric(y)
  est <- lsdv(y,x,n,t)
  est$fitted.values <- as.vector(x %*% est$coefficients)
  est$residuals <- y - est$fitted.values
  est$call <- match.call()
  class(est) <- "Lsdv"
  est
  
  
}
print.Lsdv <- function(x,...)
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
summary.Lsdv<-function(object,...)
{
  se <- sqrt(diag(object$vcov))
  tval <-coef(object)/se
  Rsqu<-object$Rsq
  TAB <- cbind(Estimate = coef(object),StdErr = se, t.value = tval, p.value = 2*pt(-abs(tval), df=object$df))
  colnames(TAB) <- c("Estimate", "Std.Err", "T value", "Pr(>z)")
  tab<-cbind(FisherStat=object$Effect,p.value=object$Peffect, R.squared=object$Rsq)
  colnames(tab) <- c("F-stat", "P-value", "R-squared")
  res <- list(call=object$call,
              coefficients=TAB, Effect=tab)
  
  res
}
#' formula
#' 
#' @param formula PIB~INF+TIR
#' @param data the dataframe 
#' @param n the number of section
#' @param t the time per section
#' @param ... not used
#' @export
Lsdv.formula <-function(formula, data=list(),n,t,...)
{
  t<-t
  n<-n
  mf <- model.frame(formula=formula, data=data)
  x <- model.matrix(attr(mf, "terms"), data=mf)
  y <- model.response(mf)
  est <- Lsdv.default(y,x,n,t,...)
  est$call <- match.call()
  est$formula <- formula
  est
}