
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

S2sls<-function(y,x,w){
  H<-cbind(x,w%*%x,w^2%*%x)
  Z<-cbind(x,w%*%y)
  zhat<-H%*%(ginv(crossprod(H))%*%crossprod(H,Z))
  beta<-ginv(crossprod(zhat,Z))%*%crossprod(zhat,y)
  yhat<-Z%*%beta
  res<-y-yhat
  resvar<-sum(res^2)/length(y)
  varcoef<-resvar*ginv(crossprod(zhat))
  std<-sqrt(diag(varcoef))
  st<-beta/std
  df<-nrow(Z)-ncol(Z)
  #pv<-2*pt(-abs(st),df=df)
  list(coefficients=beta,resvar=resvar,vcov=varcoef,std=std,st=st,df=df)
}
#' method
#'
#' @author Zaghdoudi Taha
#' @param x a numeric design matrix for the model.
#' @param ... not used
#' @export
s2sls<- function(x,...){UseMethod("s2sls") }

s2sls.default <- function(y,x,w,...)
{
  x<-as.matrix(x)
  h<-as.matrix(w)
  y<-as.numeric(y)
  est <- S2sls(y,x,w)
  #est$fitted.values <- as.vector(x %*% est$coefficients)
  est$residuals <- y - est$fitted.values
  est$call <- match.call()
  class(est) <- "s2sls"
  est


}

print.s2sls <- function(x,...)
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
#' @importFrom stats printCoefmat pt
#' @export
summary<-function(object,...)
{
  pv<-2*pt(-abs(object$st),df=object$df)
  res <- cbind(object$coefficients,object$std, object$st,pv )
  colnames(res) <- c("Estimates", "Std.Err", "T-value", "P-Value")
  cat("Formula:")
  print(object$equa)
  printCoefmat(res,has.Pvalue=TRUE)
}
#' formula
#'
#' @param formula log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp
#' @param w is the contiguity matrix
#' @param data the dataframe
#' @param ... not used
#' @importFrom stats model.frame model.matrix model.response
#' @export
s2sls.formula <-function(formula,data=list(),w,...)
{
  mf <-model.frame(formula=formula, data=data)
  x <- model.matrix(attr(mf, "terms"), data=mf)
  y <- model.response(mf)
  h<-w
  est <- s2sls.default(y,x,h,...)
  est$call <- match.call()
  est$equa <- formula
  est
}
