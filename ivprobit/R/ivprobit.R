#'Foreign-Exchange Derivatives Use By Large U.S. Bank Holding Com-panies (1996-2000).
#'
#'@docType data
#'@name eco.RData
NULL


#'Instrumental Variables Probit function
#'
#'@author Zaghdoudi Taha
#'@param  y is the dichotomous l.h.s. 
#'@param x1 is the r.h.s. exogenous variables
#'@param  y2 is the r.h.s. endogenous variables 
#'@param  x is the complete set of instruments
#'@return Coefficients: a named vector of coefficients
#'@return se: a named vector of standard error
#'@return tval: a named vector of wald statistics
#'@return pval: a named vector P-values
#'@return df: the degree of freedom
#'@export
ivprob<-function(y,x1,y2,x){
  
  y<-y
  X<-as.matrix(x)
  Y<-as.matrix(y2)
  X1<-as.matrix(x1)
  one<-matrix(1,length(y),1)
  XX<-matrix(c(one,X),ncol=ncol(X)+1)
  XX1<-matrix(c(one,X1),ncol=ncol(X1)+1)
  Z<-cbind(XX1,Y)
  d<-invpd(t(XX)%*%XX)%*%t(XX)%*%Z
  kx = ncol(X)
  ky = ncol(Y)                               
  #-----------------------------
  lin<-lm(y2~X)
  yhat<-lin$fitted.values
  uhat<-lin$residuals
  #-----------------------------
  
  d<-chol2inv(chol(t(XX)%*%XX))%*%t(XX)%*%Z
  #step2
  prob<-glm(y~X+uhat,family=binomial(link="probit"))
  J=glmvcov(prob)
  alph=prob$coefficients 
  alpha=alph[1:(kx+1)]
  lam = alph[(kx+2):(kx+ky+1)]
  Jinv=J[1:(kx+1),1:(kx+1)]
  #step3
  prob2<-glm(y~X1+uhat+yhat,family = binomial("probit"))
  bet=prob2$coefficients
  beta = bet[(length(bet)-ky+1):length(bet)]
  rho = lam-beta
  #step4
  rhoY=Y%*%rho
  ry =  rhoY
  lin2<-lm(ry~X)
  v2=vcov(lin2)
  omega =(v2+Jinv)
  #step5
  cov=invpd(t(d)%*%invpd(omega)%*%d)
  se = sqrt(diag(cov))
  delt = cov%*%t(d)%*%invpd(omega)%*%alpha
  df=nrow(x1)-ncol(x)
  #res<-cbind(delt,se,delt/se,p.value = 2*pt(-abs(delt/se), df=df))
  r1<-lm(y~X1+yhat)
  names<-rownames(coef(summary(r1)))
  #rownames(res) <- rownames(coef(summary(r1)))
  #colnames(res) <- c("Coef","S.E.","t-stat","p-val")
  #res
  list(coefficients=delt,se=se,tval=(delt/se),pval=2*pt(-abs(delt/se), df=df),df=df,names=names)
  
}

invpd<-function(mat){
  
  chol2inv(chol(mat))
  
}
glmvcov<-function(obj){
  so<-summary(obj,corr=F)
  so$dispersion * so$cov.unscaled
  
}
#' method
#' 
#' @author Zaghdoudi Taha
#' @param x a numeric design matrix for the model.
#' @param ... not used
#' @export
ivprobit <- function(x,...){UseMethod("ivprobit") }

ivprobit.default <- function(y,x1,y2,x,....){
  
  est <- ivprob(y,x1,y2,x)
  est$vcov<-est$cov
  #est$fitted.values <- as.vector(x %*% est$coefficients)
  #est$residuals <- y - est$fitted.values
  est$call <- match.call()
  class(est) <- "ivprobit"
  est
  
}
print.ivprobit <- function(x,...)
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
summary.ivprobit<-function(object,...)
{
  res<-cbind(coef(object),object$se,object$tval,object$pval)
  
  rownames(res) <- object$names
  colnames(res) <- c("Coef","S.E.","t-stat","p-val")
  printCoefmat(res,P.values=TRUE,has.Pvalue=TRUE)
  
}
#' formula
#' 
#' @param formula1 first formula
#' @param formula2 second formula
#' @param data the dataset 
#'@param ... not used
#'@export
ivprobit.formula <-function(formula1,formula2,data=list(),...)
{
  mff <- update(formula1, ~ . -1)
  mf1 <- model.frame(mff, data)
  y <- model.response(mf1)
  x1 <- model.matrix(attr(mf1, "terms"), mf1)
  # endog
  mff2 <- update(formula2, ~ . -1)
  mf2 <- model.frame(mff2, data)
  y2 <- model.response(mf2)
  x <- model.matrix(attr(mf2, "terms"), mf2)
  est <- ivprobit.default(y,x1,y2,x,...)
  est$call <- match.call()
  est$formula1 <- formula
  est
}