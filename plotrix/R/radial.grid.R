radial.grid<-function(labels=NA,label.pos=NULL,radlab=FALSE,radial.lim=NULL,
 start=0,clockwise=FALSE,label.prop=1.1,grid.pos,grid.col="gray",
 grid.bg="transparent",show.radial.grid=TRUE) {

 par(xpd=TRUE)
 if(is.null(label.pos)) label.pos<-seq(0,1.8*pi,length=9)
 if(!is.null(labels)) {
  if(is.na(labels[1])) labels<-as.character(round(label.pos,2))
 }
 if(clockwise) label.pos<--label.pos
 if(start) label.pos<-label.pos+start
 # display the circumferential grid
 angles<-seq(0,1.96*pi,by=0.04*pi)
 for(i in seq(length(grid.pos),1,by=-1)) {
  xpos<-cos(angles)*(grid.pos[i]-radial.lim[1])
  ypos<-sin(angles)*(grid.pos[i]-radial.lim[1])
  polygon(xpos,ypos,border=grid.col,col=grid.bg)
 }
 maxlength<-max(grid.pos)-radial.lim[1]
 # display the radial grid
 if(show.radial.grid) {
  xpos<-cos(label.pos)*maxlength
  ypos<-sin(label.pos)*maxlength
  segments(0,0,xpos,ypos,col=grid.col)
  xpos<-cos(label.pos)*maxlength
  ypos<-sin(label.pos)*maxlength
 }
 # display the circumferential labels
 if(!is.null(labels)) {
  xpos<-cos(label.pos)*maxlength*label.prop
  ypos<-sin(label.pos)*maxlength*label.prop
  if(radlab) {
   for(label in 1:length(labels)) {
    if(radlab < 0)
     labelsrt<-180*label.pos[label]/pi-90+
      180*(label.pos[label] > pi && label.pos[label] < 2*pi)
    else
     labelsrt<-(180*label.pos[label]/pi)+
      180*(label.pos[label] > pi/2 && label.pos[label] < 3*pi/2)
    text(xpos[label],ypos[label],labels[label],cex=par("cex.axis"),
     srt=labelsrt)
   }
  }
  else
   boxed.labels(xpos,ypos,labels,ypad=0.7,border=FALSE,
    cex=par("cex.axis"))
 }
 par(xpd=FALSE)
}
