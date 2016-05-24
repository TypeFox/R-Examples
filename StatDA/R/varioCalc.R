varioCalc <- function(X,Y,el,max.dist=300,title="",km=TRUE,plot=TRUE)
{
# Calculate variograms omnidirectional, E-W, N-S
# km ... TRUE if distances given in km, FALSE for m
# plot ... if TRUE the variogram is plotted, otherwise only the parameters are returned
#
# binned variogram
vario.b <- variog(coords=cbind(X,Y),lambda=0, data=el,max.dist=max.dist)
# computing a directional variogram
vario.0 <- variog(coords=cbind(X,Y),data=el,lambda=0,max.dist=max.dist, direction=0,tolerance=pi/8)
vario.90 <- variog(coords=cbind(X,Y),data=el,lambda=0,max.dist=max.dist, direction=pi/2,tolerance=pi/8)
ymax=max(vario.b$v,vario.0$v,vario.90$v)
if (plot){
plot(0,0,xlab="Distance [km]",ylab="Semivariogram",xlim=c(0,max.dist),
     ylim=c(0,ymax),cex.lab=1.2,type="n",xaxt="n")
if (km) axis(1,at=axTicks(1), labels=axTicks(1))
else axis(1,at=axTicks(1), labels=axTicks(1)/1000)

title(title,cex.main=1.2)
lines(vario.b,pch=1,type="p")
lines(vario.0,pch=2,type="p")
lines(vario.90,pch=3,type="p")
legend("bottomright", legend=c("Model","omnidirectional","N-S","E-W"), 
pch=c(NA,1,2,3), lty=c(1,NA,NA,NA),lwd=c(2,NA,NA,NA))
}
return(vario.b)
}

