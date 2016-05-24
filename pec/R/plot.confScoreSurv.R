##' @export 
plot.confScoreSurv <- function(x,
                               what="mean",
                               xlim,
                               ylim,
                               col,
                               lty,
                               lwd,
                               xlab="Time",
                               ylab="Confidence score",
                               legend.x="bottomright",
                               ...){
  M <- length(x$models)
  if (missing(col)) col=1:M
  if (missing(lty)) lty=1
  if (missing(lwd)) lwd=1.5
  if (length(lwd)==1) lwd <- rep(lwd,M)
  if (missing(xlim)) xlim=c(0,max(x$times))
  if (missing(ylim)) ylim=c(min(sapply(x$models,function(m)min(m$score,na.rm=TRUE))),1)
  plot(0,0,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,...)
  sumCS <- summary(x)
  switch(what,
         "mean"={
           meanScore <- sumCS[,grep("meanScore",names(sumCS))]
           loopModels <- lapply(1:M,function(m){
             lines(x$times,meanScore[,m],lty=lty[m],col=col[m],lwd=lwd[m])
           })},
         "individual"={
           loopModels <- lapply(1:M,function(m){
             loopSubjects <- apply(x$models[[m]]$score,1,function(i){
               lines(x$times,i,lty=lty[m],col=col[m],lwd=lwd[m])
             })})
           if (smooth==TRUE){
             nix <- lapply(1:length(M),function(m){
               smooth3 <- prodlim::meanNeighbors(x=x$models[[m]]$meanPred,y=x$models[[m]]$score,bandwidth=NULL)
               lines(averageY~uniqueX,data=smooth3,lty=lty[m],col=col[m],lwd=3)
             })
           }
           legend(x=legend.x,
                  xpd=NA,
                  bty="n",
                  names(x$models),
                  lty=lty,
                  col=col,
                  pch=1)
         },
         "quantiles"={
           medScore <- sumCS[,grep("medianScore",names(sumCS))]
           lowScore <- sumCS[,grep("iqrLowScore",names(sumCS))]
           upScore <- sumCS[,grep("iqrUpScore",names(sumCS))]
           loopModels <- lapply(1:M,function(m){
             lines(x$times,medScore[,m],lty=lty[m],col=col[m],lwd=lwd[m],lty=1)
             lines(x$times,lowScore[,m],lty=lty[m],col=col[m],lwd=lwd[m],lty=3)
             lines(x$times,upScore[,m],lty=lty[m],col=col[m],lwd=lwd[m],lty=3)
           })}
         )
  legend(x=legend.x,
         y=1.05,
         xpd=NA,
         bty="n",
         names(x$models),
         col=col,
         lty=lty,
         lwd=lwd[1])
}
