#' @import Rcpp
#' @import rpart.plot
#' @import survival
#' @import pec
#' @import Ecdat
#' @import rpart
#' @useDynLib DStree
#' @export



DStree<-function(formula,status,data,control=control,weights=NULL){
  Call <- match.call()
  if(is.factor(data$status)){
    stop("status variable must be numeric")
  }
  
  
  if (is.numeric(status)){
    status.col<-status
  }else
  {status.col <- which(colnames(data)==status)}  
  
  mf <- model.frame(formula,data,na.action=NULL)
  x <- model.matrix(attr(mf, "terms"), data=mf)
  y <- model.response(mf)
  dat <- data
  dat[,status.col] <- dat[,status.col] + 1
  
  rpart.resp <- paste("cbind(y,",status,")~",sep="")
  rpart.pred <- paste(attr(attributes(mf)$terms,"term.labels"),collapse="+")
  rpart.formula <- as.formula(paste(rpart.resp,rpart.pred))
  controls <- DStree.control()
  if (!missing(control)) controls[names(control)] <- control
  fit<-rpart(rpart.formula,method=alist,data=dat,control=control,
             weights=weights,y=TRUE,x=TRUE,model=TRUE)
  
  colnames(fit$cptable)[3] <- "Deviance" 
  fit$cptable <- fit$cptable[,-3]
  fit$names <- c(attributes(mf)$names[1],status)
  fit$call<-Call
  ncol.yval2<-dim(fit$frame$yval2)[2]-1
  fit$frame$yval2 <- matrix(fit$frame$yval2[,-1],ncol=ncol.yval2)
  ncol.yval2<-dim(fit$frame$yval2)[2]
  colnames(fit$frame$yval2)<-colnames(fit$frame$yval2,do.NULL=FALSE) 
  colnames(fit$frame$yval2)[1:(ncol.yval2/2)]<-paste("Haz.",1:(ncol.yval2/2),sep="")
  colnames(fit$frame$yval2)[(ncol.yval2/2+1):ncol.yval2]<-paste("Surv.",1:(ncol.yval2/2),sep="")
  class(fit)<-"DStree"
  return(fit)
}
