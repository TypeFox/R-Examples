MCtest<-function(RAN, RQ, nper=100)
{
library(gtools)
R<-list()
Q<-list()
for (i in 1:nper)
{RAN.p<-apply(RAN,2,permute)
R[[i]]<-compute.RQ(RAN.p)[,1]
Q[[i]]<-compute.RQ(RAN.p)[,2]
if (i%%100==0) cat(i,"  ")
}
R2<-as.matrix(as.data.frame(R))
colnames(R2)=NULL
Q2<-as.matrix(as.data.frame(Q))
colnames(Q2)=NULL
R.high<-rep(NA,nrow(RAN))
R.low<-rep(NA,nrow(RAN))
Q.high<-rep(NA,nrow(RAN))
Q.low<-rep(NA,nrow(RAN))

for (i in 1:nrow(RAN))
{
R.high[i]<-length(which(R2[i,]>=RQ[i,1]))/nper
R.low[i] <-length(which(R2[i,]<=RQ[i,1]))/nper
Q.high[i]<-length(which(Q2[i,]>=RQ[i,2]))/nper
Q.low[i] <-length(which(Q2[i,]<=RQ[i,2]))/nper
}
res<-cbind(R.high,R.low, Q.high, Q.low)
rownames(res)<-rownames(RAN)
return(res)
}
