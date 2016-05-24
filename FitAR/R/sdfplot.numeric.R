`sdfplot.numeric` <-
function(obj, ...){
p<-SelectModel(obj, lag.max=ceiling(length(obj)/4), Best=1)
y<-obj-mean(obj)
ans<-FitAR(y, 1:p)
PlotARSdf(ans$phiHat, InnovationVariance=ans$sigsqHat, logSdf=TRUE, ...)
}

