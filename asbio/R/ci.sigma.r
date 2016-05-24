ci.sigma<-function(data,conf=.95,S.sq=NULL,n=NULL,summarized=FALSE){
if(summarized==FALSE){
  S.sq<-var(data)
  n<-nrow(as.matrix(data))}
alpha<-1-conf
lower.chi<-qchisq(1-(alpha/2),df=n-1)
upper.chi<-qchisq((alpha/2),df=n-1)
CI<-c(S.sq,S.sq*(n-1)/lower.chi,S.sq*(n-1)/upper.chi)
head<-paste(paste(as.character(conf*100),"%",sep=""),c("Confidence interval for population variance"))
ends<-c("Estimate",paste(as.character(c((1-conf)/2,1-((1-conf)/2))*100),"%",sep=""))
res<-list(ci=CI,head=head,ends=ends)
class(res)<-"ci"
res
}
