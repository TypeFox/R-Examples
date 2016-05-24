wrapdataL <-
function(fil="sim.txt",plotte=TRUE,nInMixture=1:5,tit="Maximized at number of contributors"){
mydataL=function(x,nContributors=2) log(dataL(nContributors,p=x[,4]))
logLik=rep(NA,length(nInMixture))
dat=read.table(file=fil,header=TRUE)
dat.split=split(dat,dat[,2])
s=0
for (i in nInMixture){
s=s+1
logLikEachMarker=lapply(dat.split,mydataL,nContributors=i)
logLik[s]=sum(unlist(logLikEachMarker))
}
if (plotte){
plot(nInMixture,logLik,xlab="Num. in mixture",ylab="log(likelihood)")
title(tit)
}
oo=order(logLik,decreasing=TRUE)
res=data.frame(nInMixture[oo],logLik[oo])
colnames(res)=c("Num. in mixture","log(likelihood)")
res
}

