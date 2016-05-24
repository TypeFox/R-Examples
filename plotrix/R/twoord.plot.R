color.axis<-function(side=1,at=NULL,labels=TRUE,axlab=NA,axlab.at=NA,
 col=par("fg"),cex.axis=par("cex.axis"),cex=par("cex")) {

 xylim<-par("usr")
 if(side %% 2) {
  if(min(at) < xylim[1]) at<-at[at >= xylim[1]]
  if(max(at) > xylim[2]) at<-at[at <= xylim[2]]
 }
 else {
  if(min(at) < xylim[3]) at<-at[at >= xylim[3]]
  if(max(at) > xylim[4]) at<-at[at <= xylim[4]]
 }
 axcol<-par("col.axis")
 par(col.axis=col)
 axis(side=side,at=at,labels=rep("",length(at)),col=col)
 if(labels[1] == TRUE && length(labels) == 1) labels<-at
 mtext(labels,side=side,at=at,line=0.7,col=col,cex=cex.axis)
 par(col.axis=axcol)
 if(!is.na(axlab)) {
  if(is.na(axlab.at))
   axlab.at<-ifelse(side%%2,sum(xylim[1:2])/2,sum(xylim[3:4])/2)
  mtext(axlab,side=side,at=axlab.at,line=2,col=col,cex=cex)
 }
 if(side == 1) abline(h=xylim[3],col=col)
 if(side == 2) abline(v=xylim[1],col=col)
 if(side == 3) abline(h=xylim[4],col=col)
 if(side == 4) abline(v=xylim[2],col=col)
}

twoord.plot<-function(lx,ly,rx,ry,data=NULL,main="",xlim=NULL,lylim=NULL, 
 rylim=NULL,mar=c(5,4,4,4),lcol=1,rcol=2,xlab="",
 lytickpos=NA,ylab="",ylab.at=NA,rytickpos=NA,rylab="",rylab.at=NA,
 lpch=1,rpch=2,type="b",xtickpos=NULL,xticklab=NULL,halfwidth=0.4,
 axislab.cex=1,do.first=NULL,...) {

 if(!is.null(data)) {
  ly<-unlist(data[ly])
  ry<-unlist(data[ry])
  if(missing(lx)) lx<-1:length(ly)
  else lx<-unlist(data[lx])
  if(missing(rx)) rx <- 1:length(ry)
  else rx<-unlist(data[rx])
 }
 if(missing(lx)) lx<-1:length(ly)
 if(missing(ry)) {
  if(missing(rx)) {
   rx<-1:length(ry)
   ry<-ly
   ly<-lx
   lx<-1:length(ly)
  }
  else {
   ry<-rx
   rx<-1:length(ry)
  }
 }
 oldmar<-par("mar")
 par(mar=mar)
 if(is.null(xlim)) xlim<-range(c(lx,rx))
 if(missing(lx)) lx<-1:length(ly)
 if(is.null(lylim)) {
  lylim<-range(ly,na.rm=TRUE)
  lyspan<-diff(lylim)
  if(lyspan == 0) lyspan<-lylim[1]
  lylim[2]<-lylim[2]+lyspan*0.04
  if(lylim[1] != 0) lylim[1]<-lylim[1]-lyspan*0.04
 }
 if(length(type) < 2) type<-rep(type,2)
 # first display the "left" plot
 if(match(type[1],"bar",0)) {
  oldcex<-par(cex=axislab.cex)
  plot(lx,ly,xlim=xlim,ylim=lylim,xlab=xlab,ylab="",yaxs="i",type="n", 
   main="",axes=FALSE,...)
  par(oldcex)
  if(!is.null(do.first)) eval(parse(text=do.first))
  ybottom<-par("usr")[3]
  if (lylim[1] < 0) abline(h=0,lty=2)
  rect(lx-halfwidth,ifelse(ly<0,ly,ybottom),lx+halfwidth,
   ifelse(ly>0,ly,0),col=lcol)
 }
 else {
  oldcex<-par(cex=axislab.cex)
  plot(lx,ly,xlim=xlim,ylim=lylim,xlab=xlab,ylab="",yaxs="i",type="n", 
   main="",axes=FALSE,...)
  par(oldcex)
  if(!is.null(do.first)) eval(parse(text=do.first))
  points(lx,ly,col=lcol,pch=lpch,type=type[1],...)
 }
 title(main=main)
 xylim<-par("usr")
 #mtext(ylab,2,2,col=lcol,cex=axislab.cex)
 box()
 # display the X axis
 if(is.null(xticklab)) axis(1,cex.axis=axislab.cex)
 else {
  if(is.null(xtickpos)) xtickpos<-1:length(xticklab)
  axis(1,at=xtickpos,labels=xticklab,cex.axis=axislab.cex)
 }
 # display the left axis
 if(is.na(lytickpos[1])) lytickpos<-pretty(ly)
 if(is.na(ylab.at)) ylab.at<-mean(lytickpos)
 color.axis(2,at=lytickpos,axlab=ylab,axlab.at=ylab.at,
  col=ifelse(is.na(lcol),1,lcol),cex.axis=axislab.cex,cex=axislab.cex)
 # get the "right" y limits
 if(is.null(rylim)) {
  rylim<-range(ry,na.rm=TRUE)
  ryspan<-diff(rylim)
  if(ryspan == 0) ryspan<-rylim[1]
  rylim[2]<-rylim[2]+ryspan*0.04
  if(rylim[1] != 0) rylim[1]<-rylim[1]-ryspan*0.04
 }
 # multiplier for the "right" y values
 ymult<-diff(lylim)/diff(rylim)
# offset for the "right" y values
 yoff<-lylim[1]-rylim[1]*ymult
 if(match(type[2],"bar",0)) {
  if(rylim[1] < 0) abline("h", 0)
  rect(rx-halfwidth,ifelse(ry<0,ry,rylim[1]*ymult+yoff),rx+halfwidth,
   ifelse(ry>0,ry*ymult+yoff,0),col=rcol)
 }
 else points(rx,ry*ymult+yoff,col=rcol,pch=rpch,type=type[2],...)
 if(is.na(rytickpos[1])) rylabels<-pretty(rylim)
 else rylabels<-rytickpos
 if(min(rylabels) < rylim[1]) rylabels<-rylabels[rylabels >= rylim[1]]
 if(max(rylabels) > rylim[2]) rylabels<-rylabels[rylabels <= rylim[2]]
 axat<-rylabels*ymult+yoff
 if(is.na(rylab.at)) rylab.at<-mean(rytickpos)
 if(!is.na(rylab.at)) rylab.at<-rylab.at*ymult+yoff
 # display the right axis
 color.axis(4,at=axat,labels=rylabels,axlab=rylab,axlab.at=rylab.at,
  col=ifelse(is.na(rcol),1,rcol),cex.axis=axislab.cex,cex=axislab.cex)
 par(mar=oldmar,new=FALSE,col.axis="black")
}
