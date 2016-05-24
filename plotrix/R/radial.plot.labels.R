radial.plot.labels<-function(lengths,radial.pos=NULL,units="radians",
 radial.lim=NULL,start=0,clockwise=FALSE,labels,adj=NULL,pos=NULL,...) {

 npoints<-length(lengths)
 if(is.null(radial.pos))
  radial.pos<-seq(0,pi*(2-2/npoints),length.out=npoints)
 else {
  # clock24 starts at "midnight" and always runs clockwise
  if(units == "clock24") radial.pos<-pi*(450-radial.pos*15)/180
  # polar starts at 3 o'clock and runs counterclockwise by default
  if(units == "polar") radial.pos<-pi*radial.pos/180
 }
 if(clockwise && units != "clock24") radial.pos<--radial.pos
 if(start && units != "clock24") radial.pos<-radial.pos+start
 if(length(radial.pos) < npoints)
  radial.pos<-rep(radial.pos,length.out=npoints)
 if(is.null(radial.lim)) radial.lim<-range(lengths)
 lengths<-lengths-radial.lim[1]
 # get the vector of x positions
 xpos<-cos(radial.pos)*lengths
 # get the vector of y positions
 ypos<-sin(radial.pos)*lengths
 text(x=xpos,y=ypos,labels=labels,adj=adj,pos=pos,...)
}
