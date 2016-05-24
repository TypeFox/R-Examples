summary.nparLD<-function(object, ...)
{
 a<-object$input
 new.object<-nparLD(formula=a$formula,data=a$data,subject=a$subject,
description=FALSE,time1.order=a$time1.order,time2.order=a$time2.order,group1.order=a$group1.order, 
group2.order=a$group2.order,plot.CI=FALSE,alpha=a$alpha,show.covariance=a$show.covariance,order.warning=FALSE)

 cat("Model:","\n")
 cat(new.object$model.name,"\n","\n")
 cat("Call:","\n")
 print(a$formula)
 cat("\n","Relative Treatment Effect (RTE):","\n",sep="")
 print(new.object$RTE)
 cat("\n","Wald-Type Statistc (WTS):","\n",sep="")
 print(new.object$Wald.test)
 cat("\n","ANOVA-Type Statistc (ATS):","\n",sep="")
 print(new.object$ANOVA.test)
 if(!is.null(new.object$ANOVA.test.mod.Box))
 {
  cat("\n","Modified ANOVA-Type Statistic for the Whole-Plot Factors:","\n",sep="")
  print(new.object$ANOVA.test.mod.Box)
 }

 summarylist<-new.object
}
