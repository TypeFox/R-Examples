`sdfplot.FitAR` <-
function(obj, ...){
mti<-obj$DataTitle
sti<-obj$ModelTitle
PlotARSdf(obj$phiHat,InnovationVariance=obj$sigsqHat,logSdf=TRUE, main=mti, sub=sti, ...)
}

