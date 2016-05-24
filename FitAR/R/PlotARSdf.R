`PlotARSdf` <-
function(phi=NULL, theta=NULL, units="radial",logSdf=FALSE,InnovationVariance=1, main=NULL, sub=NULL, lwd=3, col="blue", plotQ=TRUE, ...){
sdf<-InnovationVariance
if (!is.null(phi))  
    sdf<-ARSdf(phi)*sdf
if (!is.null(theta))
    sdf<-sdf/ARSdf(theta)
if (is.null(phi)&&is.null(theta))
    sdf<-sdf*ARSdf(0)
f0<-(1/length(sdf))*(1:length(sdf))*0.5
if (units=="radial")
    f <- f0*pi
else
    f<-f0
yl<-"Spectral Density"
if (logSdf){
    sdf<-log(sdf)
    yl<-"Log Sdf"
    }
if (plotQ)    
	plot(f,sdf,type="l",xlab="frequency",ylab=yl, main=main, sub=sub, lwd=lwd, col=col, ...)
invisible(matrix(c(f,sdf),ncol=2))
}
