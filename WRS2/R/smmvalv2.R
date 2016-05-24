smmvalv2<-function(dfvec,iter=10000,alpha=.05,SEED=TRUE){
#
if(SEED)set.seed(1)
dfv<-length(dfvec)/sum(1/dfvec)
vals<-NA
tvals<-NA
J<-length(dfvec)
z=matrix(nrow=iter,ncol=J)
for(j in 1: J)z[,j]=rt(iter,dfvec[j])
vals=apply(z,1,max)
vals<-sort(vals)
ival<-round((1-alpha)*iter)
qval<-vals[ival]
qval
}