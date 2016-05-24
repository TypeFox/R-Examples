AP.test<-function(Y){
ranks<-matrix(nrow=nrow(Y),ncol=ncol(Y),data=rank(Y))
R<-apply(ranks,2,mean)
a<-ncol(Y)
b<-nrow(Y)
S<-(b-1)*cov(ranks)/(b-a+1)
C<-diag(1,nrow=a-1,ncol=a)
for(i in 1:a-1){
C[i,i+1]<--1}
F.star<-b/(a-1)*t(C%*%R)%*%(solve(C%*%S%*%t(C)))%*%C%*%R
nu1<-a-1
nu2<-(a-1)*(b-1)
AP.Table<-data.frame(df1=nu1,df2=nu2,F=F.star,p.val=pf(F.star,nu1,nu2,lower.tail=FALSE),row.names="X")
colnames(AP.Table)<-c("df1","df2","F*","P(>F)")
AP.Table
}
