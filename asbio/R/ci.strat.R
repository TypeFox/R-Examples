ci.strat<-function(data,strat,N.h,conf=.95,summarized=FALSE,use.t=FALSE,n.h=NULL,x.bar.h=NULL,var.h=NULL){
N<-sum(N.h)
k<-nrow(as.matrix(N.h))
  if(summarized==FALSE){
    var.h<-tapply(data,strat,function(x){var(x,na.rm=TRUE)})
    x.bar.h<-tapply(data,strat, function(x){mean(x,na.rm=TRUE)})
    n.h<-summary(as.factor(strat))
  }
S.xbar.str<-sqrt(sum((1-(n.h/N.h))*((N.h/N)^2)*(var.h/n.h)))
S.T.hat<-sqrt(sum((1-(n.h/N.h))*N.h^2*(var.h/n.h)))
T.hat<-sum(N.h*x.bar.h)
x.bar.str<-T.hat/N
  if(use.t==FALSE){
    z.star<-qnorm(1-((1-conf)/2))
    CI.mu<-as.matrix(c(x.bar.str,x.bar.str-z.star*S.xbar.str, x.bar.str+z.star* S.xbar.str))
    CI.T<-as.matrix(c(T.hat,T.hat-z.star* S.T.hat, T.hat+z.star* S.T.hat))
  }
  if(use.t==TRUE){
    t.star<-qt(1-((1-conf)/2),df=(sum(n.h)-k))
    CI.mu<-as.matrix(c(x.bar.str,x.bar.str-t.star*S.xbar.str, x.bar.str+t.star* S.xbar.str))
    CI.T<-as.matrix(c(T.hat,T.hat-t.star* S.T.hat, T.hat+t.star* S.T.hat))
  }
  
result<-list()
ends<-c((1-conf)/2,1-((1-conf)/2))*100
rownames(CI.mu)<-c("estimate",paste(ends[1],"%"),paste(ends[2],"%"))
rownames(CI.T)<-c("estimate",paste(ends[1],"%"),paste(ends[2],"%"))
CI<-cbind(CI.mu,CI.T)
colnames(CI)<-c("CI.mu","CI.T")
strat.summary=as.matrix(cbind(N.h,n.h,x.bar.h,var.h))
rownames(strat.summary)<-levels(strat)
result<-list(strat.summary=strat.summary,CI=t(CI))
result
}



