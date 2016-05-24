`RacfPlot` <-
function(obj, lag.max=1000, SquaredQ=FALSE, ylab=" "){
modN<-paste("Residual autocorrelation:",obj$ModelTitle,"\nSimultaneous 95% Interval")
tb<-obj$RacfMatrix
ra<-tb[,1]
sdra<-tb[,2]
MXL<-min(length(ra),lag.max)
if (SquaredQ) {
    res<-residuals(obj)
    ra<-(acf(res^2, MXL, type="correlation",plot=FALSE)$acf)[-1]
    sdra<-rep(1/sqrt(length(res)),MXL)
}
else {
    sigma<-sqrt(obj$sigsqHat)
    ra<-ra[1:MXL]
    sdra<-sdra[1:MXL]
}
clim<-qnorm(1-0.05/(2*MXL))
glim<-max(abs(ra),clim*max(sdra))*1.05
lags<-1:MXL
plot(lags, ra, type="h", ylim=c(-glim,glim), xlab="lag", ylab=ylab)
lines(c(0, length(ra)), c(0,0), col="magenta")
lines(lags, clim*sdra, col="blue")
lines(lags, -clim*sdra, col="blue")
title(main=modN) 
}

