size_n_color<-function(x=NULL,y,size,sizefun="sqrt",col,main="",
 xlim=NA,xlab="",xat=NULL,xaxlab=NULL,xcex=1,xlas=0,xgrid=FALSE,
 ylim=NA,ylab="",yat=NULL,yaxlab=NULL,ycex=1,ylas=1,ygrid=TRUE,
 mar=c(5,4,4,2),boxit=TRUE,add=FALSE,...) {

 if(!is.na(sizefun)) size<-do.call(sizefun,list(size))
 if(is.matrix(size)) {
  dimsize<-dim(size)
  if(is.null(x)) {
   x<-matrix(rep((1:dimsize[2])*max(size*2),each=dimsize[1]),ncol=dimsize[2])
   y<-matrix(rep((dimsize[1]:1)*max(size*2),dimsize[2]),ncol=dimsize[2])
   if(is.null(xat)) xat<-x[1,]
   if(is.null(yat)) yat<-y[,1]
  }
  else {
   if(is.null(xat)) xat<-1:dimsize[2]
   if(is.null(yat)) yat<-1:dimsize[1]
  }
 }
 else {
  dimsize=c(length(size),1)
  if(is.null(x)) x<-1:length(size)
 }
 xylim<-par("usr")
 aspect_ratio<-(xylim[2]-xylim[1])/(xylim[4]-xylim[3])
 maxsize<-max(size)
 if(is.na(xlim[1]))
  xlim<-c(min(x)-maxsize*aspect_ratio,max(x)+maxsize*aspect_ratio)
 if(is.na(ylim[1]))
  ylim<-c(min(y)-maxsize/aspect_ratio,max(y)+maxsize/aspect_ratio)
 if(!add) {
  oldmar<-par(mar=mar)
  plot(x,y,main=main,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,
   axes=FALSE,type="n",...)
  if(xgrid) segments(xat,xylim[3],xat,xylim[4],col="lightgray",lty=2)
  if(ygrid) segments(xylim[1],yat,xylim[2],yat,col="lightgray",lty=2)
  axis(1,at=xat,labels=xaxlab,las=xlas,cex.axis=xcex)
  axis(2,at=yat,labels=yaxlab,las=ylas,cex.axis=ycex)
  if(boxit) box()
 }
 if(is.matrix(size)) {
  if(is.null(dim(col))) col=matrix(col,nrow=dimsize[1],ncol=dimsize[2])
  for(row in 1:dimsize[1]) {
   for(column in 1:dimsize[2])
    draw.circle(x[row,column],y[row,column],size[row,column],
     border=NA,col=col[row,column])
  }
 }
 else {
  if(length(col) < length(size)) col=rep(col,size)
  for(index in 1:length(size))
   draw.circle(x[index],y[index],size[index],border=NA,col=col[index])
 }
 if(!add) par(mar=oldmar)
}
