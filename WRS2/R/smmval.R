smmval<-function(dfvec,iter=10000,alpha=.05,SEED=TRUE){
if(SEED)set.seed(1)
dfv<-length(dfvec)/sum(1/dfvec)
vals<-NA
tvals<-NA
J<-length(dfvec)
for(i in 1:iter){
for(j in 1:J){
tvals[j]<-rt(1,dfvec[j])
}
vals[i]<-max(abs(tvals))
}
vals<-sort(vals)
ival<-round((1-alpha)*iter)
qval<-vals[ival]
qval
}
