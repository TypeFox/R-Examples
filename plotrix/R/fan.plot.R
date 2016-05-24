fan.plot<-function (x,edges=200,radius=1,col=NULL,align.at=NULL, 
 max.span=NULL,labels=NULL,labelpos=NULL,label.radius=1.2,align="left", 
 shrink=0.02,main="",ticks=NULL,include.sumx=FALSE,...) {

 if(!is.numeric(x) || any(is.na(x) | x <= 0)) 
  stop("fan.plot: x values must be positive.")
 nticks<-ifelse(ticks==1,sum(x),ticks)
 if(include.sumx) x<-c(x,sum(x))
 xorder<-order(x, decreasing = TRUE)
 oldmar<-par("mar")
 # convert values into angular extents
 if (align == "center") x<-pi*x/sum(x)
 else x<-2*pi*x/sum(x)
 # stretch or shrink the values if max.span has been specified
 if(!is.null(max.span)) x<-x*max.span/max(x)
 if(is.null(align.at)) {
  if(align == "center") align.at<-pi/2
  if(align == "right") align.at<-(pi - x[xorder[1]])/2
  if(align == "left") align.at<-(pi + x[xorder[1]])/2
 }
 nx<-length(x)
 tempradius<-radius
 if(is.null(col)) col<-rainbow(nx)
 else if (length(col) < nx) col<-rep(col,length = nx)
 if(align == "left") lowpoint<-min(sin(align.at-x))
 else lowpoint<-min(sin(align.at+x))
 if(lowpoint > 0) lowpoint<-0
 par(mar=c(1,1,1,1),xpd=TRUE)
 xspan<-max(label.radius+0.012*max(nchar(labels)))
 plot(0,xlim=c(-xspan,xspan),ylim=c(lowpoint,xspan), 
  xlab="",ylab="",type="n",axes=FALSE)
 xy<-par("usr")
 pinxy<-par("pin")
 ymult<-(xy[4]-xy[3])/(xy[2]-xy[1])*(pinxy[1]/pinxy[2])
 if(nchar(main)) 
  text(0,(max(label.radius)+0.2)*ymult,main,cex=1.5,adj=c(0.5,1))
 for(i in 1:nx) {
  n<-edges*x[xorder[i]]/pi
  if(align == "center") 
   t2p<-seq(align.at-x[xorder[i]],align.at+x[xorder[i]],length=n)
  if (align == "right") 
   t2p<-seq(align.at,align.at+x[xorder[i]],length=n)
  if(align=="left") 
   t2p<-seq(align.at-x[xorder[i]],align.at,length=n)
  xc<-c(cos(t2p)*tempradius,0)
  yc<-c(sin(t2p)*tempradius*ymult,0)
  polygon(xc,yc,col=col[xorder[i]],...)
  tempradius<-tempradius-shrink
 }
 if(!is.null(ticks)) {
  if(align == "left")
   tickles<-seq(align.at,align.at-x[xorder[1]],by=-sum(x)/ticks)
  if(align == "right") 
   tickles<-seq(align.at,align.at+x[xorder[1]],by=sum(x)/ticks)
  if(align == "center") 
   tickles<-c(seq(align.at,align.at+x[xorder[1]],by=sum(x)/(2*ticks)),
    seq(align.at,align.at-x[xorder[1]],by=-sum(x)/(2*ticks)))
  stickout<-rep(c(0.04,rep(0.02,4),0.03,rep(0.02,4)),
   length.out=length(tickles))
  xpos1<-cos(tickles) * (radius + stickout)
  ypos1<-sin(tickles) * ymult * (radius + stickout)
  xpos2<-cos(tickles) * radius
  ypos2<-sin(tickles) * ymult * radius
  segments(xpos1, ypos1, xpos2, ypos2)
 }
 if(!is.null(labels)) {
  par(xpd=TRUE)
  if(align == "center") {
   labelside<-rep(c(1, -1),length.out=nx-1)
   lpos<-c(align.at+labelside*(x[xorder[1:(nx-1)]]+
    (x[xorder[2:nx]]-x[xorder[1:(nx-1)]])/2),align.at)
  }
  if(align == "right") {
   lpos<-align.at+(x[xorder]-c((x[xorder[1:(nx-1)]]-
    x[xorder[2:nx]])/2,x[xorder[nx]]/2))
   labelside<-rep(1,nx)
  }
  if(align == "left") {
   lpos<-align.at-(x[xorder]-c((x[xorder[1:(nx-1)]]-
    x[xorder[2:nx]])/2,x[xorder[nx]]/2))
   labelside<-rep(-1,nx)
  }
  if(is.null(labelpos)) {
   labelpos<-lpos
   ldiff<-abs(diff(labelpos))
   squeezers<-which(ldiff<0.15)
   if(length(squeezers)) {
    for(squeeze in squeezers) {
     labelpos[1:squeeze]<-labelpos[1:squeeze]+
      (0.15-ldiff[squeeze]) * labelside[1:squeeze]
     labelpos[(squeeze+1):nx]<-labelpos[(squeeze+1):nx]-
      (0.15-ldiff[squeeze]) * labelside[(squeeze+1):nx]
    }
   }
  }
  innerend<-seq(0.98,by=-shrink,length=nx)
  endpos<-lpos
  if(length(label.radius) < nx) 
   label.radius<-rep(label.radius,length.out=nx)
  xpos1<-cos(labelpos)*label.radius[xorder]
  ypos1<-sin(labelpos)*ymult*label.radius[xorder]
  for(i in 1:nx)
   boxed.labels(xpos1[i],ypos1[i],labels[xorder[i]], 
    border=FALSE,xpad=0,ypad=0.1,adj=c(0.5-cos(labelpos[i])/2.5,0.5))
  xpos1<-cos(labelpos)*(label.radius[xorder]+radius)/2
  ypos1<-sin(labelpos)*ymult*(label.radius[xorder]+radius)/2
  xpos2<-cos(endpos)*radius*innerend
  ypos2<-sin(endpos)*ymult*radius*innerend
  segments(xpos1,ypos1,xpos2,ypos2)
 }
 par(mar=oldmar,xpd=FALSE)
 invisible(labelpos)
}
