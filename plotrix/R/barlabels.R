barlabels<-function(xpos,ypos,labels=NULL,cex=1,prop=0.5,miny=0,offset=0,...) {
 if(is.data.frame(ypos)) ypos<-as.matrix(ypos)
 if(is.null(labels)) labels<-ypos
 # usually don't want to display zero labels
 display<-ypos > miny
 if(is.matrix(ypos)) {
  # prop is within the scope of the current environment
  cumcenter<-function(x,pos) return(cumsum(x)-x*prop)
  stacked<-length(xpos) < length(ypos)
  if(stacked) {
   # replicate the x positions one by one, but the offsets group by group
   xpos<-rep(xpos,each=length(ypos)/length(xpos))+
    rep(c(-offset,offset),length(ypos)/(2*length(xpos)))
   ypos<-apply(ypos,2,cumcenter)
  }
  else ypos<-ypos*prop
 }
 else ypos<-ypos*prop
 # allow labels to extend beyond the plot area
 par(xpd=TRUE)
 boxed.labels(xpos[display],ypos[display],labels[display],cex=cex,...)
 par(xpd=FALSE)
}
