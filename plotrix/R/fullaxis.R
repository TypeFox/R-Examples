fullaxis<-function(side=1,at=NULL,labels=TRUE,line=NA,pos=NA,outer=FALSE,
 font=NA,lty="solid",lwd=1,lwd.ticks=lwd,col=NULL,col.ticks=NULL,
 hadj=NA,padj=NA,...) {

 xylim<-par("usr")
 if(!is.na(pos)) xylim[c(3,1,4,2)[side]]<-pos
 par(xpd=TRUE)
 switch(side,
  segments(xylim[1],xylim[3],xylim[2],xylim[3],col=col,lty=lty,lwd=lwd),
  segments(xylim[1],xylim[3],xylim[1],xylim[4],col=col,lty=lty,lwd=lwd),
  segments(xylim[1],xylim[4],xylim[2],xylim[4],col=col,lty=lty,lwd=lwd),
  segments(xylim[2],xylim[3],xylim[2],xylim[4],col=col,lty=lty,lwd=lwd))
 par(xpd=FALSE)
 axpos<-axis(side=side,at=at,labels=labels,line=line,pos=pos,outer=outer,
  font=font,lty=lty,lwd=lwd,col=col,col.ticks=col.ticks,
  hadj=hadj,padj=padj,...)
 invisible(axpos)
}
