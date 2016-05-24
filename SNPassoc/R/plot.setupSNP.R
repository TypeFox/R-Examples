`plot.setupSNP`<-
function(x,which=1,...){
   plot(x[,attr(x,"colSNPs")[which]], label=labels(x)[which], ...)
}