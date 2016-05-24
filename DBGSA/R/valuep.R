valuep<-function(MM,N,pname){
library(fdrtool)
N<-as.matrix(N);
n=nrow(MM);
m=ncol(N);
E<-matrix(nrow=nrow(MM),ncol=3);
for(i in 1:n)
{
w=0;
  for(j in 1:m){
               if(MM[i,3]>N[MM[i,2],j])
                  {
                   w=w+1;
                   }
                     }
p=1-w/m
E[i,1]=p
}

fdr=fdrtool(E[,1],statistic="pvalue")
E[,2]<-fdr$qval
E[,3]<-fdr$lfdr
T<-data.frame(GOid=MM[,1],GeneNum=MM[,2],pvalue=E[,1],Fdr=E[,2],local_fdr=E[,3])
write.csv(T,file=pname,row.names=FALSE)
}