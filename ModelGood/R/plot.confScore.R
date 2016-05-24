plot.confScore <- function(x,
                           original=TRUE,
                           smooth=FALSE,
                           xlim,
                           ylim,
                           col,
                           lwd,
                           xlab,
                           ylab="Confidence score",
                           legend.x="bottomleft",
                           ...){
  if (missing(col)) col=1:length(x$models)
  if (missing(xlim)) xlim=c(0,1)
  if (missing(ylim)) ylim=c(min(sapply(x$models,function(m)min(m$score,na.rm=TRUE))),1)
  if (missing(xlab)){
    if (original)
      xlab="Model prediction"
    else
      xlab="Mean bootstrap prediction"}
  plot(0,0,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,type="n",...)
  nix <- lapply(1:length(x$models),function(m){
    if (original)
      f <- score~origPred
    else
      f <- score~meanPred
    points(f,data=x$models[[m]],col=col[m])
  })
  if (smooth==TRUE){
    nix <- lapply(1:length(x$models),function(m){
      smooth3 <- prodlim::meanNeighbors(x=x$models[[m]]$meanPred,y=x$models[[m]]$score,bandwidth=NULL)
      lines(averageY~uniqueX,data=smooth3,col=col[m],lwd=3)
    })
  }
  legend(x=legend.x,xpd=NA,bty="n",names(x$models),col=col,pch=1)
}
