plot.nparLD<-function(x, ...)
{
 a<-x$input
 new.object<-nparLD(formula=a$formula,data=a$data,subject=a$subject,
description=a$description,time1.order=a$time1.order,time2.order=a$time2.order,group1.order=a$group1.order, 
group2.order=a$group2.order,plot.CI=TRUE,alpha=a$alpha,show.covariance=a$show.covariance,order.warning=FALSE)

 summarylist<-new.object
}
