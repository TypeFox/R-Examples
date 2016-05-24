`sdfplot.Arima` <-
function(obj, ...){
sti<-paste("ARMA(",obj$arma[1],",",obj$arma[2],")",sep="")
if (obj$arma[5]>1){
    sti<-paste("S",sti, "(",obj$arma[3],",",obj$arma[4],")",obj$arma[5],sep="")
}
PlotARSdf(phi=obj$model$phi, theta=-(obj$model$theta), InnovationVariance=obj$sigma2, logSdf=TRUE, sub=sti, ...)
}

