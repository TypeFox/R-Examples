
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

threesls<-function(y,reg,iv){

  #================ 2sls==============
  N<-length(y)
  df<-nrow(reg)-ncol(reg)
  th<-crossprod(iv,iv)
  thx<-crossprod(iv,reg)
  xhat<-iv%*%ginv(th)%*%thx
  beta<-ginv(crossprod(xhat,reg))%*%crossprod(xhat,y)
  yhat<-reg%*%beta
  res<-y-yhat
  resvar<-sum(res^2)/(N)
  varcoef<-resvar*ginv(crossprod(xhat,xhat))
  sd<-sqrt(diag(varcoef))
  st<-beta/sd
  rsq<-(cor(yhat,y))^2
  #==================3sls=============
  sig<-crossprod(res)
  sig<-sig/(length(y))
  bb<-kronecker(sig,iv)
  b<-crossprod(iv,bb)
  b<-ginv(b)
  a<-crossprod(iv,reg)
  c<-crossprod(iv,y)
  d<-ginv(t(a)%*%b%*%a)
  beta<-d%*%crossprod(a,b)%*%c
  pval<- 2*pt(-abs(st),df=df)
  list(beta=beta,sd=sd,st=st,pval=pval,rsq=rsq)
}
#' method
#'
#' @author Zaghdoudi Taha
#' @param x a numeric design matrix for the model.
#' @param ... not used
#' @export
tsls <- function(x,...){UseMethod("tsls") }



tsls.default <- function(y,x,h,...)
{


  x<-as.matrix(x)
  h<-as.matrix(h)
  y<-as.numeric(y)
  est <- threesls(y,x,h)
  #est$fitted.values <- as.vector(x %*% est$coefficients)
  #est$residuals <- y - est$fitted.values
  est$call <- match.call()
  class(est) <- "tsls"
  est


}
print.tsls <- function(x,...)
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
summary.tsls<-function(object,...)
{
  result <- cbind(object$beta,object$sd, object$st,object$pval)
  colnames(result) <- c("Estimates", "Std.Err", "T-value", "P-Value")
  cat("Formula:")
  print(object$equa)
  cat("\nRsquared :",object$rsq,"\n")
  printCoefmat(result,has.Pvalue=TRUE)
}
#' formula
#'
#' @param formula PIB~INF+TIR|Cap+m2r
#' @param data the dataframe
#' @param ... not used
#' @import Formula
#' @export
tsls.formula <-function(formula,data=list(),...)
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
  est <- tsls.default(y,x,h,...)
  est$call <- match.call()
  est$equa <- formula
  est
}
