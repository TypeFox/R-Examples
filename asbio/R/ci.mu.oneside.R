ci.mu.oneside<-function(data,conf=.95,n=NULL,Var=NULL,xbar=NULL,summarized=FALSE,N=NULL,fpc = FALSE,tail="upper"){
  if(summarized==TRUE)Var.xbar<-Var/n
  if(summarized==FALSE){
    xbar<-mean(data)
    n<-nrow(as.matrix(data))
    Var.xbar<-var(data)/n}
S.xbar<-ifelse(fpc==FALSE,sqrt(Var.xbar),sqrt((1-(n/N))*Var.xbar))
tstar<-qt(conf,df=n-1)
  if(tail=="upper")CI<-c(xbar,NA,xbar+tstar*S.xbar)
  if(tail=="lower")CI<-c(xbar,xbar-tstar*S.xbar,NA)
res<-list()
res$head<-paste(paste(as.character(conf*100),"%",sep=""),c("One-sided confidence interval for population mean"))
res$ci<-CI
res$ends<-c("Estimate",paste(as.character(c((1-conf),1-((1-conf)))*100),"%",sep=""))
class(res)<-"ci"
res
}

