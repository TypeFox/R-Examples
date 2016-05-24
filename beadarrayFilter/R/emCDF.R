emCDF <-
function(iccResults,iccQuant=1)
{
cols<-brewer.pal(12,"Paired")
if(length(iccQuant)==1)
{
plot(ecdf(iccResults$icc[,2]),ylim=c(0.3,1),xlim=c(0,1),ylab="Fraction of non-informative bead types",xlab=expression(paste("Intracluster correlation threshold (",rho,")")),main="",col=cols[1])
legend("bottomright", legend=c(paste("q", iccQuant,  sep = "")) ,col=cols[1],lty=1)
lines(c(0.5,0.5),c(0,1),lty=2) 
text(0.5,0.95,expression(threshold: rho>=0.5),pos=4)
}
else 
{
allicc<-reshape(iccResults$icc, idvar="ProbeID", varying=list(names(iccResults$icc[,-1])),v.names="ICC", direction="long")
plot(ecdf(allicc$ICC[allicc$time==1]),ylim=c(0.3,1),xlim=c(0,1),ylab="Fraction of non-informative bead types",xlab=expression(paste("Intracluster correlation threshold (",rho,")")),main="")
for (i in unique(allicc$time))
{
lines(ecdf(allicc$ICC[allicc$time==i]),col=cols[i])
legend("bottomright", legend=c(paste("q", iccQuant,  sep = "")) ,col=cols[1:length(unique(allicc$time))],lty=1) 
lines(c(0.5,0.5),c(0,1),lty=2) 
text(0.5,0.95,expression(threshold: rho>=0.5),pos=4)
}
}
}
