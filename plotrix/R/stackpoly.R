stackpoly<-function(x,y=NULL,main="",xlab="",ylab="",xat=NA,xaxlab=NA,
 xlim=NA,ylim=NA,lty=1,lwd=1,border=NA,col=NULL,staxx=FALSE,stack=FALSE,
 axis2=TRUE,axis4=TRUE,padj=0,...) {

 ydim<-dim(y)
 if(is.null(y[1])) {
  y<-x
  ydim<-dim(y)
  if(is.null(ydim)) x<-1:length(y)
  else x<-matrix(rep(1:ydim[1],ydim[2]),ncol=ydim[2])
 }
 if(stack) y<-t(unlist(apply(as.matrix(y),1,cumsum)))
 if(is.na(xlim[1])) xlim<-range(x)
 if(is.na(ylim[1])) ylim<-range(y)
 plot(0,main=main,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,type="n",
  xaxs="i",yaxs="i",axes=FALSE,...)
 if(is.matrix(y) || is.list(y)) {
  plotlim<-par("usr")
  if(is.na(xat[1])) xat<-x[,1]
  if(is.na(xaxlab[1])) xaxlab<-xat
  if(staxx) staxlab(at=xat,labels=xaxlab)
  else axis(1,at=xat,labels=xaxlab,padj=padj)
  if(axis2) axis(2)
  if(axis4) axis(4)
  if(is.null(col[1])) col=rainbow(ydim[2])
  else if(length(col)<ydim[2]) col<-rep(col,length.out=ydim[2])
  if(length(border) < ydim[2]) border<-rep(border,length.out=ydim[2])
  if(length(lty)<ydim[2]) lty<-rep(lty,length.out=ydim[2])
  if(length(lwd)<ydim[2]) lwd<-rep(lwd,length.out=ydim[2])
  for(pline in seq(ydim[2],1,by=-1)) {
   if(pline==1) {
    polygon(c(x[1],x[,pline],x[ydim[1]]),
     c(plotlim[3],y[,pline],plotlim[3]),
     border=border[pline],col=col[pline],lty=lty[pline],lwd=lwd[pline])
   }
   else
    polygon(c(x[,pline],rev(x[,pline-1])),c(y[,pline],rev(y[,pline-1])),
     border=border[pline],col=col[pline],lty=lty[pline],lwd=lwd[pline])
  }
 }
 else {
  polygon(c(min(x),x,max(x),0),c(0,y,0,0),border=border,col=col,lty=lty,
   lwd=lwd)
  if(is.na(xat[1])) xat<-x
  if(is.na(xaxlab[1])) xaxlab<-xat
  if(staxx) staxlab(at=xat,labels=xaxlab)
  else axis(1,at=xat,labels=xaxlab)
  if(axis2) axis(2)
  if(axis4) axis(4)
 }
 box()
}
