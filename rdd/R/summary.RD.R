#' Summarizing Regression Discontinuity Designs
#' 
#' \code{summary} method for class \code{"RD"}
#' 
#' @method summary RD
#' @param object an object of class \code{"RD"}, usually a result of a call to \code{\link{RDestimate}}
#' @param digits number of digits to display
#' @param ... unused
#' @return \code{summary.RD} returns an object of \link{class} "\code{summary.RD}" which has the following components:
#' \item{coefficients}{A matrix containing bandwidths, number of observations, estimates, 
#' SEs, z-values and p-values for each estimated bandwidth.}
#' \item{fstat}{A global F-test of the corresponding model}
#' @include RDestimate.R
#' @importFrom stats residuals pf
#' @export
#' @author Drew Dimmery <\email{drewd@@nyu.edu}>

summary.RD<-function(object,digits=max(3, getOption("digits") - 3),...){
  cat("\n")
  cat("Call:\n")
  print(object$call)
  cat("\n")
  
  cat("Type:\n")
  cat(object$type,"\n\n")
  
  #If the model wasn't included in the output, we need to get it
  mod<-FALSE
  if("model" %in% names(object$call)) mod<-object$call$model
  if(!mod){
    object$call$model<-TRUE
    object$call$verbose<-FALSE
    object<-eval.parent(object$call)
  }
  n<-length(object$est)

  obs<-vector(length=n)
  if(object$type=="sharp") {
    for(i in 1:n) obs[i]<-length(residuals(object$model[[i]]))
  } else {
    for(i in 1:n) obs[i]<-length(residuals(object$model$iv[[i]]))
  }
  cat("Estimates:\n")
  #Need to get this to give at least as much as stata does in fuzzy designs
  stars<-vector(length=n)

  for(i in 1:n) {
    stars[i]<-if (object$p[i]<0.001) "***" else if(object$p[i]<0.01) "**" else if(object$p[i]<0.05) "*" else if(object$p[i]<0.1) "." else " "
  }

  out<-matrix(
    c(object$bw,
      object$obs,
      object$est,
      object$se,
      object$z,
      object$p),
    nrow=n)
  rownames(out)<-names(object$est)
  colnames(out)<-c("Bandwidth","Observations","Estimate","Std. Error","z value","Pr(>|z|)")

  print.default(cbind(apply(out,2,function(x) format(x,digits=digits)),
                      " "=stars),quote=FALSE,print.gap=2,right=FALSE)
  cat("---\n")
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
  
  out<-list(coefficients=out)
  class(out)<-"summary.RD"
  fstat<-matrix(NA,nrow=n,ncol=4)
  if(object$type=="sharp") {
    for(i in 1:n) {
      fstat[i,]<-c(summary(object$model[[i]])$fstatistic[[1]],
               summary(object$model[[i]])$fstatistic[[2]],
               summary(object$model[[i]])$fstatistic[[3]],
        pf(summary(object$model[[i]])$fstatistic[[1]],
          summary(object$model[[i]])$fstatistic[[2]],
          summary(object$model[[i]])$fstatistic[[3]])
      )
    }
    fstat[,4]<-2*apply(cbind(fstat[,4],1-fstat[,4]),1,min)
  } else {
    for(i in 1:n) {
      fstat[i,]<-c(
        summary(object$model$iv[[i]])$waldtest[1],
        summary(object$model$iv[[i]])$waldtest[3],
        summary(object$model$iv[[i]])$waldtest[4],
        summary(object$model$iv[[i]])$waldtest[2])
    }
  }
  colnames(fstat)<-c("F","Num. DoF","Denom. DoF","p")
  rownames(fstat)<-names(object$est)
  cat("F-statistics:\n")
  print.default(apply(fstat,2,function(x) format(x,digits=digits)),quote=FALSE,print.gap=2,right=FALSE)
  out$fstat<-fstat
  return(invisible(out))
}
