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
#' Instrumental least square dummy variable
#' 
#' @author Zaghdoudi Taha
#' @param y vector of dependent variable
#' @param x matrix of independents variables
#' @return Coefficients a named vector of coefficients
#' @return vcov covariance matrix of coefficients
#' @return df the degree of freedom
#' @examples
#' pib<-as.matrix(c(12,3,4,0.4,0.7,5,0.7,0.3,0.6,89,7,8,45,7,4,5,0.5,5),nrows=18,ncols=1)
#' tir<-as.matrix(c(12,0.3,4,0.4,7,12,3.0,6.0,45,7.0,0.8,44,65,23,4,6,76,9),nrows=18,ncols=1)
#' inf<-as.matrix(c(1.2,3.6,44,1.4,0.78,54,0.34,0.66,12,0.7,8.0,12,65,43,5,76,65,8),nrows=18,ncols=1)
#' npl<-as.matrix(c(0.2,3.8,14,2.4,1.7,43,0.2,0.5,23,7.8,88,36,65,3,44,65,7,34),nrows=18,ncols=1)
#' # create data frame
#' mdata<-data.frame(p=pib,ti=tir,int=inf,np=npl)
#' # create the designed matrix for the model
#' d<-matrix(c(mdata$p,mdata$int,mdata$np),18, 3)
#' dd<-matrix(c(mdata$ti,mdata$int,mdata$np),18, 3)
#' # Fit an Instrumental least square dummy variable model
#' ls<-lsdv(mdata$p,d,dd)
#' ls
tosls<-function(y,x,h){
  x<-x
  y<-y
  h<-h
  N<-length(y)
  df<-nrow(x)-ncol(x)
  th<-crossprod(h,h)
  thx<-crossprod(h,x)
  xhat<-h%*%ginv(th)%*%thx
  beta<-ginv(crossprod(xhat,x))%*%crossprod(xhat,y)
  yhat<-x%*%beta
  res<-y-yhat
  resvar<-sum(res^2)/(N)
  varcoef<-resvar*ginv(crossprod(xhat,xhat))
  scr<-sum(res^2)
  sct<-crossprod(y,y)
  rsq<-1-scr/sct
  list(coefficients=beta,vcov=varcoef,Rsq=rsq, df=df)
}

#' method
#' 
#' @author Zaghdoudi Taha
#' @param x a numeric design matrix for the model.
#' @param ... not used
#' @export
Tosls <- function(x,...){UseMethod("Tosls") }



Tosls.default <- function(y,x,h,...)
{
  
  
  x<-as.matrix(x)
  h<-as.matrix(h)
  y<-as.numeric(y)
  est <- tosls(y,x,h)
  est$fitted.values <- as.vector(x %*% est$coefficients)
  est$residuals <- y - est$fitted.values
  est$call <- match.call()
  class(est) <- "Tosls"
  est
  
  
}
print.Tosls <- function(x,...)
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
summary<-function(object,...)
{
  se <- sqrt(diag(object$vcov))
  tval <-coef(object)/se
  Rsqu<-object$Rsq
  df<-object$df
  TAB <- cbind(Estimate = coef(object),StdErr = se, t.value = tval, p.value = 2*pt(-abs(tval),df=df))
  colnames(TAB) <- c("Estimate", "Std.Err", "T value", "Pr(>z)")
  
  res <- list(call=object$call,
              coefficients=TAB)
  
  res
}
#' formula
#' 
#' @param formula PIB~INF+TIR|Cap+m2r
#' @param data the dataframe 
#' @param ... not used
#' @import Formula
#' @export
Tosls.formula <-function(formula,data=list(),...)
{
  
  
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
  est <- Tosls.default(y,x,h,...)
  est$call <- match.call()
  est$formula <- formula
  est
}