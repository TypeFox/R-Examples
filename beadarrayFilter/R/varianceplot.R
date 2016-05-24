varianceplot <-
function(iccResults,q=1,delta=0.5)
{
icc<-iccResults$icc
##there is only 1 ID variable and it is the first column
informativeData<-subset(icc, icc[,q+1]>=delta)[,1]
vars<-iccResults$withinvar
allvar<-reshape(vars, idvar="ProbeID", varying=list(names(vars[,-1])),v.names="withinvar", direction="long")
allvar$informID<-ifelse(allvar[,1] %in% informativeData,1,0)
taus<-iccResults$betweenvar
varaall<-data.frame(allvar,taus[,2])
colnames(varaall)<-c(colnames(allvar),"betweenvar")

####plot for informative max ICC filtering
plot(varaall$betweenvar,varaall$withinvar,xlim = c(0, max(varaall$betweenvar)), ylim = c(0,max(varaall$withinvar)), 
       xlab = "variance of random intercept",ylab="within array variance",type="n",cex=1.3)
points(varaall$betweenvar[varaall$informID==0],varaall$withinvar[varaall$informID==0], pch=21,col="blue")
points(varaall$betweenvar[varaall$informID==1],varaall$withinvar[varaall$informID==1], pch=24,col="red")
#legend(2, 12, c("noninformative", "informative") ,col=c("blue","red"), pch=c(21,24)) 
legend((min(varaall$withinvar) + 2), (max(varaall$withinvar) - 4), c("noninformative", "informative") ,col=c("blue","red"), pch=c(21,24)) 
}
