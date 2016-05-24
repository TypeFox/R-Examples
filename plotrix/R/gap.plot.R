# Try to rewrite this for an arbitrary number of gaps

gap.plot<-function(x,y,gap,gap.axis="y",bgcol="white",breakcol="black",
 brw=0.02,xlim=range(x),ylim=range(y),xticlab,xtics=NA,yticlab,ytics=NA,
 lty=rep(1,length(x)),col=rep(par("col"),length(x)),pch=rep(1,length(x)),
 add=FALSE,stax=FALSE,...) {

 if(missing(y) && !missing(x)) {
  y<-x
  x<-1:length(y)
 }
 if(missing(gap)) stop("gap must be specified")
 gapsize<-diff(gap)
 xaxl<-par("xlog")
 yaxl<-par("ylog")
 if(gap.axis == "y") {
  if(length(gap) > 3) ylim[2]<-ylim[2] - (gapsize[1] + gapsize[3])
  else ylim[2]<-ylim[2]-gapsize[1]
 }
 if(gap.axis == "x") {
  if(length(gap) > 3) xlim[2]<-xlim[2] - (gapsize[1] + gapsize[3])
  else xlim[2]<-xlim[2]-gapsize[1]
 }
 rangexy <- c(range(xlim),range(ylim))
 xgw<-(rangexy[2]-(rangexy[1]+gapsize))*brw
 ygw<-(rangexy[4]-(rangexy[3]+gapsize))*brw
 if(is.na(xtics[1])) xtics<-pretty(x)
 if(is.na(ytics[1])) ytics<-pretty(y)
 if(missing(xticlab)) xticlab<-xtics
 if(missing(yticlab)) yticlab<-ytics
 if(length(col) < length(y)) col<-rep(col,length.out=length(y))
 if(gap.axis == "y") {
  littleones<-which(y < gap[1])
  if(length(gap) > 3) {
   middleones<-which(y >= gap[2] + ygw & y < gap[3])
   bigones<-which(y >= gap[4] + ygw)
   lostones<-sum(c(y > gap[1] & y < gap[2] + ygw,y > gap[3] & y < gap[4] + ygw))
  }
  else {
   middleones<-NA
   bigones<-which(y >= gap[2] + ygw)
   lostones<-sum(y > gap[1] & y < gap[2] + ygw)
  }
  if(lostones) warning("some values of y may not be displayed")
 }
 else {
  littleones<-which(x < gap[1])
  if(length(gap) > 3) {
   middleones<-which(x >= gap[2] + xgw & x < gap[3])
   bigones<-which(x >= gap[4] + xgw)
   lostones<-sum(c(x > gap[1] & x < gap[2] + xgw,x > gap[3] & x < gap[4] + xgw))
   if(missing(xlim)) xlim<-c(min(x),max(x) - (gapsize[1] + gapsize[3]))
  }
  else {
   middleones<-NA
   bigones<-which(x >= gap[2])
   lostones<-sum(x > gap[1] & x < gap[2] + xgw)
   if(missing(xlim)) xlim<-c(min(x),max(x) - gapsize[1])
  }
  if(lostones) warning("some values of x will not be displayed")
  if(missing(ylim)) ylim<-range(y)
 }
 if(length(lty) < length(x)) lty<-rep(lty,length.out=length(x))
 if(length(col) < length(x)) col<-rep(col,length.out=length(x))
 if(length(pch) < length(x)) pch<-rep(pch,length.out=length(x))
 if(add) {
  points(x[littleones],y[littleones],lty=lty[littleones],
   col=col[littleones],pch=pch[littleones],...)
  if(gap.axis == "y") {
   if(length(gapsize) > 2) {
    points(x[middleones],y[middleones]-gapsize[1],
     lty=lty[middleones],col=col[middleones],pch=pch[middleones],...)
    points(x[bigones],y[bigones] - (gapsize[1] + gapsize[3]),
     lty=lty[bigones],col=col[bigones],pch=pch[bigones],...)
   }
   else points(x[bigones],y[bigones]-gapsize[1],
    lty=lty[bigones],col=col[bigones],pch=pch[bigones],...)
  }
  else {
   if(length(gapsize) > 2) {
    points(x[middleones] - gapsize[1],y[middleones],
     lty=lty[middleones],col=col[middleones],pch=pch[middleones],...)
    points(x[bigones] - (gapsize[1] + gapsize[3]),y[bigones],
     lty=lty[bigones],col=col[bigones],pch=pch[bigones],...)
   }
   else points(x[bigones]-gapsize[1],y[bigones],
    lty=lty[bigones],col=col[bigones],pch=pch[bigones],...)
  }
 }
 else {
  plot(x[littleones],y[littleones],xlim=xlim,ylim=ylim,axes=FALSE,
   lty=lty[littleones],col=col[littleones],pch=pch[littleones],...)
  box()
  if(gap.axis == "y") {
   if(!is.na(xtics[1])) axis(1,at=xtics,labels=xticlab)
   littletics<-which(ytics < gap[1])
   if(length(gapsize) > 2) {
    middletics<-which(ytics >= gap[2]+ygw & ytics <= gap[3])
    bigtics<-which(ytics >= gap[4]+ygw)
    show.at<-c(ytics[littletics],ytics[middletics] - gapsize[1],
     ytics[bigtics]-(gapsize[1] + gapsize[3]))
    show.labels<-c(yticlab[littletics],yticlab[middletics],yticlab[bigtics])
   }
   else {
    bigtics<-which(ytics >= gap[2])
    show.at<-c(ytics[littletics],ytics[bigtics] - gapsize[1])
    show.labels<-c(ytics[littletics],yticlab[bigtics])
   }
   if(stax) {
    axis(2,at=show.at,labels=rep("",length(show.labels)))
    staxlab(2,at=show.at,labels=show.labels)
   }
   else axis(2,at=show.at,labels=show.labels)
   axis.break(2,gap[1],style="gap",bgcol=bgcol,
    breakcol=breakcol,brw=brw)
   if(length(gapsize) > 2) {
    axis.break(2,gap[3]-gapsize[1],style="gap",bgcol=bgcol,
     breakcol=breakcol,brw=brw)
    points(x[middleones],y[middleones]-gapsize[1],
     lty=lty[middleones],col=col[middleones],pch=pch[middleones],...)
    points(x[bigones],y[bigones]-(gapsize[1]+gapsize[3]),
     lty=lty[bigones],col=col[bigones],pch=pch[bigones],...)
   }
   else points(x[bigones],y[bigones]-gapsize[1],
    lty=lty[bigones],col=col[bigones],pch=pch[bigones],...)
  }
  # x gaps need to be fixed
  else {
   if(!is.na(ytics[1])) axis(2,at=ytics,labels=yticlab)
   littletics<-which(xtics < gap[1])
   if(length(gapsize) > 2) {
    middletics<-which(xtics >= gap[2] + xgw & xtics <= gap[3])
    bigtics<-which(xtics > gap[4]+xgw)
    show.at<-c(xtics[littletics],xtics[middletics]-gapsize[1],
     xtics[bigtics]-(gapsize[1]+gapsize[3]))
    show.labels<-c(xticlab[littletics],xticlab[middletics],xticlab[bigtics])
   }
   else {
    bigtics<-which(xtics > gap[2]+xgw)
    show.at<-c(xtics[littletics],xtics[bigtics]-gapsize[1])
    show.labels<-c(xticlab[littletics],xticlab[bigtics])
   }
   if(stax) {
    axis(1,at=show.at,labels=rep("",length(show.labels)))
    staxlab(1,at=show.at,labels=show.labels)
   }
   else axis(1,at=show.at,labels=show.labels)
   axis.break(1,gap[1],style="gap")
   if(length(gapsize) > 2) {
    axis.break(1,gap[3]-gapsize[1],style="gap")
    points(x[middleones]-gapsize[1],y[middleones],
     lty=lty[middleones],col=col[middleones],pch=pch[middleones],...)
    points(x[bigones]-(gapsize[1]+gapsize[3]),y[bigones],
     lty=lty[bigones],col=col[bigones],pch=pch[bigones],...)
   }
   else points(x[bigones]-gapsize[1],y[bigones],
    lty=lty[bigones],col=col[bigones],pch=pch[bigones],...)
  }
 }
}
