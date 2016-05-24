get.triprop<-function(use.percentages=FALSE,cnames=c("1st","2nd","3rd")) {
 cat("Enter the label and ")
 cat(ifelse(use.percentages,"percentages ","proportions "))
 cat("of",cnames[1],cnames[2],"and",cnames[3],"for each observation.\n")
 cat("Enter a blank observation label to end.\n")
 nextlabel<-"dummy"
 nextprop<-0
 proplabels<-NA
 prop1<-NA
 prop2<-NA
 prop3<-NA
 nprop<-0
 totprop<-ifelse(use.percentages,100,1)
 tolerance<-ifelse(use.percentages,1,0.01)
 while(nchar(nextlabel)) {
  nextlabel<-readline("Observation label - ")
  if(nchar(nextlabel)) {
   if(is.na(proplabels[1])) proplabels<-nextlabel
   else proplabels<-c(proplabels,nextlabel)
   cat(cnames[1],"- ")
   nextprop<-as.numeric(readline())
   if(is.na(prop1[1])) prop1<-nextprop
   else prop1<-c(prop1,nextprop)
   cat(cnames[2],"- ")
   nextprop<-as.numeric(readline())
   if(is.na(prop2[1])) prop2<-nextprop
   else prop2<-c(prop2,nextprop)
   cat(cnames[3],"- ")
   nextprop<-as.numeric(readline())
   if(is.na(prop3[1])) prop3<-nextprop
   else prop3<-c(prop3,nextprop)
   nprop<-nprop+1
  }
  sumprop<-prop1[nprop]+prop2[nprop]+prop3[nprop]
  if(abs(totprop-sumprop) > tolerance)
   cat("Warning - sum not equal to",totprop,"\n")
 }
 triprop<-cbind(prop1,prop2,prop3)
 rownames(triprop)<-proplabels
 colnames(triprop)<-cnames
 return(triprop)
}

triax.abline<-function(b=NULL,r=NULL,l=NULL,col=par("col"),lty=par("lty"),
 cc.axes=FALSE) {

 sin60<-sin(pi/3)
 if(!is.null(b)) {
  if(any(b>1)) b<-b/100
  if(cc.axes) {
   bx2<-b+0.5*(1-b)
   by2<-sin60*(1-b)
   segments(b,0,bx2,by2,col=col,lty=lty)
  }
  else {
   bx2<-0.5*(1-b)
   by2<-sin60*(1-b)
   segments(1-b,0,bx2,by2,col=col,lty=lty)
  }
 }
 if(!is.null(r)) {
  if(any(r>1)) r<-r/100
  if(cc.axes) {
   rx1<-0.5*r
   ry1<-sin60*r
   rx2<-1-r*0.5
   segments(rx1,ry1,rx2,ry1,col=col,lty=lty)
  }
  else {
   rx1<-0.5*(r+1)
   ry1<-sin60*(1-r)
   rx2<-1-r
   segments(rx1,ry1,r,0,col=col,lty=lty)
  }
 }
 if(!is.null(l)) {
  if(any(l>1)) l<-l/100
  if(cc.axes) {
   lx1<-0.5-l*0.5
   lx2<-1-l
   ly<-sin60*(1-l)
   segments(lx1,ly,lx2,0,col=col,lty=lty)
  }
  else {
   lx1<-l*0.5
   ly<-l*sin60
   segments(lx1,ly,1-lx1,ly,col=col,lty=lty)
  }
 }
}

triax.points<-function(x,show.legend=FALSE,label.points=FALSE,
 point.labels=NULL,col.symbols=par("fg"),pch=par("pch"),
 bg.symbols=par("bg"),cc.axes=FALSE,...) {

 if(dev.cur() == 1)
  stop("Cannot add points unless the triax.frame has been drawn")
 if(missing(x))
  stop("Usage: triax.points(x,...)\n\twhere x is a 3 column array of proportions or percentages")
 if(!is.matrix(x) && !is.data.frame(x))
  stop("x must be a matrix or data frame with at least 3 columns and one row.")
 if(any(x > 1) || any(x < 0)) {
  if(any(x < 0))
   stop("All proportions must be between zero and one.")
  if(any(x > 100))
   stop("All percentages must be between zero and 100.")
  # convert percentages to proportions
  x<-x/100
 }
 if(any(abs(rowSums(x)-1) > 0.01))
  warning("At least one set of proportions does not equal one.")
 sin60<-sin(pi/3)
 if(cc.axes) {
  ypos<-x[,3]*sin60
  xpos<-x[,1]+x[,3]*0.5
 }
 else {
  ypos<-x[,3]*sin60
  xpos<-1-(x[,1]+x[,3]*0.5)
 }
 nobs<-dim(x)[1]
 points(x=xpos,y=ypos,pch=pch,col=col.symbols,bg=bg.symbols,type="p",...)
 if(is.null(point.labels)) point.labels<-rownames(x)
 if(label.points) thigmophobe.labels(xpos,ypos,point.labels)
 if(show.legend) {
  legend(0.2,0.7,legend=point.labels,pch=pch,col=col.symbols,
   xjust=1,yjust=0)
 }
 invisible(list(x=xpos,y=ypos))
}

triax.frame<-function(at=seq(0.1,0.9,by=0.1),axis.labels=NULL,
 tick.labels=NULL,col.axis="black",cex.axis=1,cex.ticks=1,
 align.labels=TRUE,show.grid=FALSE,col.grid="gray",lty.grid=par("lty"),
 cc.axes=FALSE) {

  sin60<-sin(pi/3)
  # bottom ticks
  bx1<-at
  bx2<-bx1+0.01-0.02*cc.axes
  by1<-rep(0,9)
  if(cc.axes) by2<-rep(-0.02*sin(2*pi/3),9)
  else by2<-rep(-0.02*sin60,9)
  # left ticks
  ly1<-at*sin60
  lx1<-bx1*0.5
  lx2<-lx1-0.02+0.013*cc.axes
  if(cc.axes) ly2<-ly1+rep(0.014*sin60,9)
  else ly2<-ly1
  # right ticks
  rx1<-at*0.5+0.5
  rx2<-rx1+0.01
  ry1<-rev(ly1)
  if(cc.axes) ry2<-ry1
  else ry2<-rev(ly2)+0.02*sin60
  if(show.grid) {
   par(fg=col.grid)
   segments(bx1,by1,lx1,ly1,lty=lty.grid)
   segments(lx1,ly1,rev(rx1),rev(ry1),lty=lty.grid)
   segments(rx1,ry1,bx1,by1,lty=lty.grid)
  }
  par(fg=col.axis,xpd=TRUE)
  if(is.null(tick.labels)) {
   if(cc.axes) tick.labels<-list(l=rev(at),r=rev(at),b=rev(at))
   else tick.labels<-list(l=at,r=at,b=at)
  }
  else {
   if(cc.axes) {
    tick.labels$l<-rev(tick.labels$l)
    tick.labels$r<-rev(tick.labels$r)
    tick.labels$b<-rev(tick.labels$b)
   }
  }
  # left axis label
  if(align.labels) par(srt=60)
  text(0.13,0.5,axis.labels[3-cc.axes],adj=0.5,cex=cex.axis)
  # left axis tick labels
  if(cc.axes) {
   par(srt=300)
   xoffset<-0.02
   yoffset<-0.04
  }
  else {
   par(srt=0)
   xoffset<-0.05
   yoffset<-0
  }
  text(lx1-xoffset,ly1+yoffset,tick.labels$l,cex=cex.ticks)
  # right axis label
  if(align.labels) {
   par(srt=300)
   label.adj<-0.5
  }
  else {
   par(srt=0)
   label.adj<-0
  }
  text(0.86,0.52,axis.labels[2+cc.axes],adj=label.adj,cex=cex.axis)
  # right axis tick labels
  if(cc.axes) {
   par(srt=0)
   xoffset<-0.033
   yoffset<-0.005
  }
  else {
   par(srt=60)
   xoffset<-0.015
   yoffset<-0.045
  }
  text(rx2+xoffset,ry1+yoffset,tick.labels$r,cex=cex.ticks)
  # bottom axis tick labels
  if(cc.axes) {
   par(srt=60)
   xoffset<- -0.03
  }
  else {
   par(srt=300)
   xoffset<-0.03
  }
  text(bx1+xoffset,by1-0.05,rev(tick.labels$b),cex=cex.ticks)
  # bottom axis label
  par(srt=0)
  text(0.5,-0.14,axis.labels[1],cex=cex.axis)
  # draw the triangle and ticks
  x1<-c(0,0,0.5)
  x2<-c(1,0.5,1)
  y1<-c(0,0,sin60)
  y2<-c(0,sin60,0)
  par(fg=col.axis)
  segments(x1,y1,x2,y2)
  # bottom ticks
  segments(bx1,by1,bx2,by2)
  # left ticks
  segments(lx1,ly1,lx2,ly2)
  # right ticks
  segments(rx1,ry1,rx2,ry2)
}

triax.fill<-function(col) {
 nrows<-length(col)
 sin60<-sin(pi/3)
 rowlen<-1
 xinc<-0.5/nrows
 yinc<-sin60/nrows
 for(trirow in 1:nrows) {
  startx<-0.5-xinc*(trirow-1)
  starty<-sin60-trirow*yinc
  dir<-1
  for(triangle in 1:rowlen) {
   polygon(c(startx-xinc,startx,startx+xinc),
    c(starty,starty+yinc*dir,starty),border=NA,
    col=col[[trirow]][triangle])
   startx<-startx+xinc
   starty<-starty+yinc*dir
   dir<--dir
  }
  rowlen<-rowlen+2
 }
}

triax.plot<-function (x=NULL,main="",at=seq(0.1,0.9,by=0.1),
  axis.labels=NULL,tick.labels=NULL,col.axis="black",
  cex.axis=1,cex.ticks=1,align.labels=TRUE,show.grid=FALSE,
  col.grid="gray",lty.grid=par("lty"),cc.axes=FALSE,
  show.legend=FALSE,label.points=FALSE,point.labels=NULL,
  col.symbols="black",pch=par("pch"),mar=c(5,2,4,2),no.add=TRUE,...) {

  oldpar<-par("fg","pty","mar","srt","xpd")
  par(xpd=TRUE)
  if(is.null(axis.labels)) axis.labels<-colnames(x)[1:3]
  par(pty="s",mar=mar)
  plot(0.5,type="n",axes=FALSE,xlim=c(0,1),ylim=c(0,1),main=main,
   xlab="",ylab="")
  triax.frame(at=at,axis.labels=axis.labels,
   tick.labels=tick.labels,col.axis=col.axis,cex.axis=cex.axis,
   cex.ticks=cex.ticks,align.labels=align.labels,show.grid=show.grid,
   col.grid=col.grid,lty.grid=lty.grid,cc.axes=cc.axes)
  if(is.null(x)) xypos <- NULL
  else
   xypos<-triax.points(x,show.legend=show.legend,
    label.points=label.points,point.labels=point.labels,
    col.symbols=col.symbols,pch=pch,cc.axes=cc.axes,...)
  if(no.add) par(oldpar)
  invisible(list(xypos=xypos,oldpar=oldpar))
}
