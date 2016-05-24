# Arguments:
# u,v - the x (longitude) and y (latitude) offsets
# OR orientation and magnitude
# xpos,ypos - the centers of the vectors
# scale - the longest arrow as a proportion of the cell size
# headspan - the extent of the arrowhead as a proportion of cell size
# this function doesn't assume a 1:1 aspect ratio

vectorField<-function(u,v,xpos=NA,ypos=NA,scale=1,headspan=0.1,
 vecspec=c("lonlat","rad","deg"),col=par("fg")) {
 udim<-dim(u)
 if(is.na(xpos[1])) xpos<-col(u)
 if(is.na(ypos[1])) ypos<-udim[1]-row(u)+1
 ymult<-getYmult()
 # if long/lat not specified
 if(match(vecspec[1],"lonlat",0) == 0) {
  # convert the degrees to radians if necessary
  if(match(vecspec[1],"deg",0)) u<-pi*u/180.
  # save the magnitudes
  mag<-v
  # get the x offsets without losing the direction
  tempu<-v*cos(u)
  v<-v*sin(u)*ymult
  u<-tempu
 }
 else mag<-sqrt(u*u+v*v)
 # make sure that date/time x values are their numeric equivalents
 if(is.null(dim(xpos))) maxmag<-0.5*max(diff(as.numeric(xpos)))/max(mag)
 else maxmag<-0.5*max(diff(as.numeric(xpos[1,])))/max(mag)
 u2<-u*scale*maxmag
 v2<-v*scale*maxmag
 if(is.null(udim)) length=headspan
 else length<-headspan*par("pin")[1]/udim[2]
 arrows(xpos-u2,ypos-v2,xpos+u2,ypos+v2,length=length,col=col)
}
