superpc.plotred.lrtest<- function(object.lrtestred, call.win.metafile=FALSE){
  
if(call.win.metafile) {win.metafile()}

n.components=object.lrtestred$n.components

na.range=function(x){range(x[!is.na(x)])}

  if(!is.null(object.lrtestred$lrtest.reduced.lower)){
ylim=na.range(c(object.lrtestred$lrtest.reduced.lower, object.lrtestred$lrtest.reduced.upper))
}
else{ylim=na.range(object.lrtestred$lrtest.reduced)}

if(n.components==2){par(mfrow=c(1,2))}
if(n.components==3){par(mfrow=c(2,2))}
par(mar=c(6,4,5,2))

for(j in 1:n.components){
ylab=paste("LR statistic, #Components=",as.character(j),sep="")
plot(object.lrtestred$shrinkages, object.lrtestred$lrtest.reduced[,j],xlab="Shrinkage amount",ylab=ylab, ylim=ylim,type="b")

abline(h = qchisq(0.95, j), lty = 2)

if(!is.null(object.lrtestred$lrtest.reduced.lower)){
error.bars(object.lrtestred$shrinkages, object.lrtestred$lrtest.reduced.lower[,j],
object.lrtestred$lrtest.reduced.upper[,j], lty=2)
}
axis(3,at=object.lrtestred$shrinkages, labels=as.character(object.lrtestred$num.features[,j]), cex=0.7)
mtext("Number of genes", 3, 4, cex = 0.8)
}


  

if(call.win.metafile){dev.off()}
return(TRUE)

}

error.bars <-function(x, upper, lower, width = 0.005, ...) {
  xlim <- range(x)
  barw <- diff(xlim) * width
  segments(x, upper, x, lower, ...)
  segments(x - barw, upper, x + barw, upper, ...)
  segments(x - barw, lower, x + barw, lower, ...)
  range(upper, lower)
}


