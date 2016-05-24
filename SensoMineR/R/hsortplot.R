hsortplot=function(don,group,numr=2,numc=2){
nb.juge=length(group)
gp=0
for(i in 1:nb.juge){
distance=group[i]-tab.disjonctif(don[,(gp+1):(gp+group[i])])%*%t(tab.disjonctif(don[,(gp+1):(gp+group[i])]))
cah=cluster::agnes(distance,method="single",diss=TRUE)
if (((i-1)%%(numc*numr))==0){
dev.new()
par(mfrow = c(numr, numc))
}
plot(cah,main=paste("H-Sorting of subject",i,sep=" "),xlab="",sub="",which.plots=2)
gp=gp+group[i]
}}