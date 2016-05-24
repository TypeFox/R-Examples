battleship.plot<-function(x,mar=c(2,5,5,1),col="white",border="black", 
 main="",xlab="",ylab="",xaxlab=NA,yaxlab=NA,cex.labels=1,
 maxxspan=0.45,maxyspan=0.45) {

 dimx<-dim(x)
 if(length(dimx) != 2) 
  stop("battleship.plot(x) where x is a 2 dimensional matrix or data frame")
 if(is.data.frame(x)) x<-as.matrix(x)
 oldmar<-par("mar")
 par(mar=mar)
 plot(0,xlim=c(0.5,dimx[2]+0.5),ylim=c(0.5,dimx[1]+0.5),axes=FALSE,
  type="n",xlab="",ylab="")
 title(main=main,line=mar[3]-2)
 mtext(xlab,side=1,line=0)
 if(is.na(xaxlab[1])) {
  xaxlab<-colnames(x)
  if(is.null(xaxlab)) xaxlab<-1:dimx[2]
 }
 staxlab(side=3,at=1:dimx[2],labels=xaxlab,srt=45,adj=0,top.line=0,
  ticklen=0,cex=cex.labels)
 if(is.na(yaxlab[1])) {
  yaxlab<-rownames(x)
  if(is.null(yaxlab)) yaxlab<-1:dimx[1]
 }
 staxlab(side=2,at=dimx[1]:1,labels=yaxlab,nlines=1,srt=0,adj=1,
  top.line=0,ticklen=0,cex=cex.labels)
 normx<-maxxspan*x/max(x,na.rm=TRUE)
 rect(rep(1:dimx[2],each=dimx[1])-normx,rep(dimx[1]:1,dimx[2])-maxyspan,
  rep(1:dimx[2],each=dimx[1])+normx,rep(dimx[1]:1,dimx[2])+maxyspan,
  col=col,border=border)
 par(mar=oldmar)
}
