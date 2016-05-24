plot.summary.rmas.risk <-
function(x, ylim=NULL, col=NULL,xlab=NULL,
                                   ylab=NULL, main=NULL,...){
  cosa <- x
  if(is.null(main)) main<-cosa$main   
  if(is.null(ylab)) ylab<-"probability"
  if(is.null(xlab)) xlab<-"population threshold"
  if(is.null(ylim)) ylim<- range(cosa[,-1])
  if(is.null(col)) col<- c("blue", "red", "red")
  plot(cosa$Threshold, cosa$Probability, type="l", ylim=ylim, col=col[1],
       ylab=ylab, xlab=xlab, main=main,...)
  lines(cosa$Threshold, cosa[,3], col=col[2], lty=2)
  lines(cosa$Threshold, cosa[,4], col=col[3], lty=2)
}

