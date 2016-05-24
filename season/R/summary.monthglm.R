summary.monthglm<-function(object, ...){
  if (class(object)!="monthglm"){stop("Object must be of class 'monthglm'")} 
## Tabulate the monthly data ##
  z<-qnorm(0.975)
  s<-summary(object$glm)
  type<-as.character(object$call$family)[1]
  out<-as.data.frame(matrix(data=NA,nrow=nrow(s$coef),ncol=5))
  names(out)<-c('mean','lower','upper','zvalue','pvalue')
  row.names(out)<-row.names(s$coef)
  out$mean<-s$coef[,1]
  out$lower<-s$coef[,1]-(z*s$coef[,2])
  out$upper<-s$coef[,1]+(z*s$coef[,2])
  out$zvalue<-s$coef[,3]
  out$pvalue<-s$coef[,4]
# Exponentiate the results if rate or odds ratio
  if (type=="poisson"|type=="binomial"){
     out$mean<-exp(out$mean)
     out$lower<-exp(out$lower)
     out$upper<-exp(out$upper)
  }
# Just keep results with months
  index<-grep("months",row.names(out),ignore.case=TRUE,value=FALSE)
  totable<-out[index,] # Select months
  effect<-''
  if (type=="poisson"){effect='RR'}
  if (type=="binomial"){effect='OR'}
# returns
 ret<-list()
 ret$n=length(object$residuals)
 ret$month.ests=totable
 ret$month.effect=effect
 class(ret) <- "summary.monthglm"
 ret # uses print.summary.monthglm
}
