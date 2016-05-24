feather.plot<-function(r,theta,xpos,yref=0,use.arrows=TRUE,
 col.refline="lightgray",fp.type="s",main="",xlab="",ylab="",
 xlabels=NULL,...) {

 if(missing(xpos)) xpos<-1:length(theta)
 # reverse the angles and start at 12 o'clock
 if(fp.type == "m") theta<-5*pi/2-theta
 # calculate the x offset of the end of the lines
 x<-r*cos(theta)
 # and the y offset
 y<-r*sin(theta)
 # get the aspect ratio of the plot
 xmult<-diff(range(xpos))/(diff(range(y))*2)
 # correct the x offsets
 #x<-x*xmult
 xlim<-range(c(xpos,x+xpos))
 ylim<-range(c(y,yref))
 xdiff<-xlim[2]-xlim[1]
 ydiff<-ylim[2]-ylim[1]
 # equalize the ranges of x and y
 if(xdiff > ydiff) {
  addon<-(xdiff-ydiff)/2
  ylim<-ylim+c(-addon,addon)
 }
 if(ydiff > xdiff) {
  addon<-(ydiff-xdiff)/2
  xlim<-xlim+c(-addon,addon)
 }
 oldpty<-par(pty="s")
 plot(0,xlim=xlim,ylim=ylim,type="n",main=main,xlab=xlab,ylab=ylab,
  axes=FALSE)
 box()
 if(is.null(xlabels)) axis(1)
 else axis(1,at=xpos,labels=xlabels)
 abline(h=yref,col=col.refline)
 if(use.arrows) arrows(xpos,yref,xpos+x,y+yref,length=0.1,...)
 else segments(xpos,yref,xpos+x,y+yref,...)
 par(pty=oldpty)
}
