barp<-function(height,width=0.4,names.arg=NULL,legend.lab=NULL,legend.pos=NULL,
 col=NULL,border=par("fg"),main=NULL,xlab="",ylab="",xlim=NULL,ylim=NULL,x=NULL,
 staxx=FALSE,staxy=FALSE,height.at=NULL,height.lab=NULL,cex.axis=par("cex.axis"),
 pch=NULL,cylindrical=FALSE,shadow=FALSE,do.first=NULL,ylog=FALSE,srt=NA) {

 height.class<-attr(height,"class")
 if(!is.null(height.class)) {
  if(match(height.class,"dstat",0)) {
   md1<-length(height)
   md2<-dim(height[[1]])[2]
   meanmat<-matrix(NA,nrow=md1,ncol=md2)
   colnames(meanmat)<-colnames(height[[1]])
   for(row in 1:md1) meanmat[row,]<-height[[row]][1,]
   height<-meanmat
  }
  if(match(height.class,"freq",0)) height<-height[[1]] 
 }
 if(is.data.frame(height)) its_ok<-is.numeric(unlist(height))
 else its_ok<-is.numeric(height)
 if(!its_ok) stop("barp can only display bars with numeric heights")
 hdim<-dim(height)
 if(is.null(x)) x<-1:length(height)
 if(is.null(hdim) || length(hdim) == 1) {
  ngroups<-length(height)
  barcol<-col
  barpinfo<-list(x=x,y=height)
  hdim<-NULL
 }
 else {
  ngroups<-hdim[2]
  x<-1:ngroups
  if(!is.matrix(col) && length(col)==hdim[1])
   barcol<-matrix(rep(col,each=ngroups),nrow=hdim[1],byrow=TRUE)
  else barcol<-col
  if(!is.matrix(pch) && length(pch)==hdim[1])
   pch<-matrix(rep(pch,ngroups),nrow=hdim[1])
  barpinfo<-list(x=matrix(rep(1:ngroups,each=hdim[1]),ncol=hdim[2]),
   y=as.matrix(height))
 }
 if(is.null(xlim)) xlim<-range(x)+c(-0.6,0.6)
 if(is.null(ylim)) {
  negy<-any(height<0,na.rm=TRUE)
  if(negy) miny<-min(height,na.rm=TRUE)*1.05
  else miny<-ifelse(ylog,min(height)/10,0)
  ylim<-c(miny,max(height,na.rm=TRUE)*1.05)
 }
 else {
  miny<-ylim[1]
  negy<-miny<0
 }
 plot(ylim[1],type="n",main=main,xlab=xlab,ylab=ylab,axes=FALSE,xlim=xlim,
  ylim=ylim,xaxs="i",yaxs="i",log=ifelse(ylog,"y",""))
 if(!is.null(do.first)) eval(parse(text=do.first))
 if(negy) abline(h=0)
 if(is.null(names.arg)) names.arg<-x
 if(staxx) {
  axis(1,at=x,labels=rep("",ngroups),cex.axis=cex.axis)
  staxlab(1,at=x,labels=names.arg,cex=cex.axis,srt=srt)
 }
 else axis(1,at=x,labels=names.arg,cex.axis=cex.axis)
 if(is.null(height.at)) {
  if(ylog) height.at<-axTicks(2,log=TRUE)
  else height.at<-pretty(ylim)
  if(max(height.at) > max(height,na.rm=TRUE))
   height.at<-height.at[-length(height.at)]
 }
 if(is.null(height.lab)) height.lab<-height.at
 if(staxy) {
  axis(2,at=height.at,labels=rep("",length(height.lab)),cex.axis=cex.axis)
  staxlab(2,at=height.at,labels=height.lab,cex=cex.axis,srt=srt)
 }
 else axis(2,at=height.at,labels=height.lab,cex.axis=cex.axis)
 bottoms<-ifelse(negy,0,miny)
 if(is.null(hdim)) {
  if(shadow) {
   for(bar in 1:ngroups)
    polygon.shadow(c(x[bar]-width,x[bar]-width,x[bar]+width,x[bar]+width),
     c(bottoms,height[bar],height[bar],bottoms),
     offset=c(0.2*width,0.05*(height[bar]-ylim[2])))
  }
  if(cylindrical)
   cylindrect(x-width,bottoms,x+width,height,col=barcol,
    border=border)
  else {
   if(is.null(pch))
    rect(x-width,bottoms,x+width,height,col=barcol,border=border)
   else
    rectFill(x-width,bottoms,x+width,height,bg="white",fg="black",
     pch=pch)
  }
 }
 else {
  bottoms<-matrix(bottoms,nrow=hdim[1],ncol=hdim[2])
  barwidth<-2*width/hdim[1]
  for(subgroup in 1:hdim[1]) {
   barpinfo$x[subgroup,]<-1:ngroups-width+(subgroup-0.5)*barwidth
   if(shadow) {
    for(bar in 1:ngroups) {
     barleft<-bar-width+(subgroup-1)*2*width/hdim[1]
     barright<-barleft+2*width/hdim[1]
     polygon.shadow(c(barleft,barleft,barright,barright),
      c(bottoms[bar],height[subgroup,bar],height[subgroup,bar],bottoms[bar]),
      offset=c(0.2*width,0.05*(height[subgroup,bar]-ylim[2])))
    }
   }
   if(cylindrical)
    cylindrect(1:ngroups-width+(subgroup-1)*barwidth,bottoms[subgroup,],
     1:ngroups-width+(subgroup)*barwidth,height[subgroup,],
     col=barcol[subgroup,],border=border)
   else {
    if(is.null(pch))
     rect(1:ngroups-width+(subgroup-1)*barwidth,bottoms[subgroup,],
      1:ngroups-width+(subgroup)*barwidth,height[subgroup,],
      col=barcol[subgroup,],border=border)
    else
     rectFill(1:ngroups-width+(subgroup-1)*barwidth,bottoms[subgroup,],
      1:ngroups-width+(subgroup)*barwidth,height[subgroup,],
      bg="white",fg="black",pch=pch[subgroup,])
   }
  }
 }
 if(!is.null(legend.lab)) {
  xjust<-yjust<-0.5
  if(is.null(legend.pos)) {
   cat("Click at the lower left corner of the legend\n")
   legend.pos<-locator(1)
   xjust<-yjust<-0
  }
  legend(legend.pos,legend=legend.lab,fill=col,xjust=xjust,yjust=yjust)
 }
 box()
 invisible(barpinfo)
}
