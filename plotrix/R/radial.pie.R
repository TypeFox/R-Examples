drawSectorAnnulus<-function(angle1,angle2,radius1,radius2,col,angleinc=0.03) {
 if(angle1 > angle2) {
  temp<-angle1
  angle1<-angle2
  angle2<-temp
 }
 if(radius1 > radius2) {
  temp<-radius1
  radius1<-radius2
  radius2<-temp
 }
 angles<-seq(angle1,angle2,by=angleinc)
 angles[length(angles)]<-angle2
 xpos<-c(cos(angles)*radius1,cos(rev(angles))*radius2)
 ypos<-c(sin(angles)*radius1,sin(rev(angles))*radius2)
 polygon(xpos,ypos,col=col,border=col)
}

# plots sectors composed of one or more sectors of annuli on a circular grid.
# radial.extents are the radial extents of the sectors (a vector),
# optionally with sub-extents that define the annuli
# (a matrix with each sector a column)
# sector.edges are the positions of the radii that define the sectors,
# defaulting to n equal sectors filling the circle where
# n is the number of radial extents.
# sector.values are the values in each sector of an annulus that
# will be represented by colors and must be the same data type and
# dimension as radial.extents.
# If a list of sector colors is passed, it will take precedence
# and sector.colors will not be scaled from sector.values

radial.pie<-function(radial.extents,sector.edges=NULL,
 sector.colors=NULL,cs1=c(0,1),cs2=c(0,1),cs3=c(0,1),alpha=1,
 labels=NA,label.pos=NULL,radlab=FALSE,start=0,clockwise=FALSE,
 label.prop=1.1,radial.lim=NULL,main="",xlab="",ylab="",mar=c(2,2,3,2),
 show.grid=TRUE,show.grid.labels=4,show.radial.grid=TRUE,
 grid.col="gray",grid.bg="transparent",grid.unit=NULL,
 radial.labels=NULL,boxed.radial=TRUE,add=FALSE,...) {
 
 if(is.null(radial.lim)) radial.lim<-range(radial.extents)
 if(is.null(sector.edges)) {
  if(clockwise)
   sector.edges<-seq(2*pi+start,start,length.out=length(radial.extents)+1)
  else
   sector.edges<-seq(start,2*pi+start,length.out=length(radial.extents)+1)
 }
 if(is.null(label.pos))
  label.pos<-sector.edges[-length(sector.edges)]+diff(sector.edges)/2
 if(show.grid) {
  if(length(radial.lim) < 3) grid.pos<-pretty(radial.lim)
  else grid.pos<-radial.lim
  if(grid.pos[1] < radial.lim[1]) grid.pos<-grid.pos[-1]
  maxlength<-max(grid.pos-radial.lim[1])
 }
 else {
  grid.pos<-NA
  maxlength<-diff(radial.lim)
 }
 oldpar<-par("xpd","mar","pty")
 if(!add) {
  par(mar=mar,pty="s")
  maxrad<-max(unlist(radial.extents))
  plot(0,xlim=c(-maxrad,maxrad),ylim=c(-maxrad,maxrad),type="n",axes=FALSE)
  if(show.grid)
   radial.grid(labels=labels,label.pos=label.pos,radlab=radlab,
    radial.lim=radial.lim,start=start,clockwise=clockwise,
    label.prop=label.prop,grid.pos=grid.pos,
    grid.col=grid.col,grid.bg=grid.bg)
 }
 fadeColor<-function(col,nfades) {
  rgbcol<-col2rgb(col)
  redinc<-(255-rgbcol[1])/nfades
  reds<-(rgbcol[1]+0:nfades*redinc)/255
  greeninc<-(255-rgbcol[2])/nfades
  greens<-(rgbcol[2]+0:nfades*greeninc)/255
  blueinc<-(255-rgbcol[3])/nfades
  blues<-(rgbcol[3]+0:nfades*blueinc)/255
  return(rgb(reds[1:nfades],greens[1:nfades],blues[1:nfades]))
 }
 nsectors<-length(radial.extents)
 if(is.list(radial.extents)) {
  if(is.null(sector.colors)) sector.colors<-rainbow(nsectors)
  for(sector in 1:nsectors) {
   annuli<-radial.extents[[sector]]
   annulus.colors<-fadeColor(sector.colors[[sector]],length(annuli))
   for(annulus in 1:(length(annuli)-1)) {
    drawSectorAnnulus(sector.edges[[sector]],sector.edges[[sector+1]],
     annuli[annulus],annuli[annulus+1],annulus.colors[annulus])    
   }
  }
 }
 else {
  if(is.null(sector.colors)) sector.colors<-rainbow(nsectors)
  for(sector in 1:nsectors) {
   drawSectorAnnulus(sector.edges[sector],sector.edges[sector+1],
    0,radial.extents[sector],sector.colors[sector])
  }
 }
 if(show.grid.labels) {
  if(show.grid.labels%%2) {
   ypos<-grid.pos-radial.lim[1]
   xpos<-rep(0,length(grid.pos))
   if(show.grid.labels==1) ypos<--ypos
  }
  else {
   xpos<-grid.pos-radial.lim[1]
   ypos<-rep(0,length(grid.pos))
   if(show.grid.labels==2) xpos<--xpos
  }
  if(is.null(radial.labels)) radial.labels<-grid.pos
  if(!is.null(grid.unit))
   radial.labels[length(grid.pos)]<-
    paste(radial.labels[length(grid.pos)],grid.unit)
  if(boxed.radial)
   boxed.labels(xpos,ypos,radial.labels,border=FALSE,
    cex=par("cex.lab"))
  else text(xpos,ypos,radial.labels,cex=par("cex.lab"))
 }
 invisible(oldpar)
}
