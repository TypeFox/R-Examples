index.AMMI<-function(model)
{
A<-model$biplot[,1:4]
A<-A[A[,1]=="GEN",-c(1,2)]
pc<-model$analysis[1,4]/model$analysis[2,4]
ASV<-apply(A,1,function(x) sqrt(pc*(x[1])^2+(x[2])^2))
rk<-rank(ASV)
B<-model$means
W<-tapply.stat(B[,3],B[,2],function(x) mean(x,rm.na=TRUE))
Rx<-rank(-W[,2])
YSI<-rk+Rx
ranking<-data.frame(ASV,YSI,rASV=rk,rYSI=Rx,means=W[,2])
invisible(ranking)
}
