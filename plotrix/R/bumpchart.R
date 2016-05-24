bumpchart<-function(y,top.labels=colnames(y),labels=rownames(y),
 rank=TRUE,mar=c(2,8,5,8),pch=19,col=par("fg"),lty=1,lwd=1,
 arrows=FALSE,...) {

 if(missing(y)) stop("Usage: bumpchart(y,top.labels,labels,...)")
 ydim<-dim(y)
 if(is.null(ydim)) stop("y must be a matrix or data frame")
 oldmar<-par("mar")
 par(mar=mar)
 if(rank) y<-apply(y,2,rank)
 # to get things the right way round, reverse the order of everything
 labels<-rev(labels)
 pch=rev(pch)
 col=rev(col)
 lty=rev(lty)
 lwd=rev(lwd)
 y<-apply(y,2,rev)
 if(arrows) {
  matplot(t(y),ylab="",type="p",pch=pch,col=col,axes=FALSE,...)
  for(row in 1:(ydim[2]-1))
   p2p_arrows(rep(row,ydim[1]),y[,row],rep(row+1,ydim[1]),y[,row+1],
    col=col,lty=lty,lwd=lwd,...)
 }
 else 
  matplot(t(y),ylab="",type="b",pch=pch,col=col,lty=lty,
   lwd=lwd,axes=FALSE,...)
 par(xpd=TRUE)
 xylim<-par("usr")
 minspacing<-strheight("M")*1.5
 text(1:ydim[2],xylim[4],top.labels)
 labelpos<-spreadout(y[,1],minspacing)
 text(xylim[1],labelpos,labels,adj=1)
 labelpos<-spreadout(y[,ydim[2]],minspacing)
 text(xylim[2],labelpos,labels,adj=0)
 par(mar=oldmar,xpd=FALSE)
}
