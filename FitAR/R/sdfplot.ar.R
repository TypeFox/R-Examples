`sdfplot.ar` <-
function(obj, ...){
sti<-paste("AR(",obj$order,")",sep="")
PlotARSdf(obj$ar, InnovationVariance=obj$var.pred, logSdf=TRUE, sub=sti, ...)
}

