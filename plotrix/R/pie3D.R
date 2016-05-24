pie3D.labels<-function(radialpos,radius=1,height=0.1,theta=pi/6, 
 labels,labelcol=par("fg"),labelcex=1.5,labelrad=1.25,minsep=0.3){

 oldcex<-par("cex")
 nlab<-length(labels)
 par(cex=labelcex,xpd=TRUE)
 for (i in 1:nlab) {
  if(i < nlab) {
   labelsep<-radialpos[i+1] - radialpos[i]
   if(labelsep < minsep) {
    radialpos[i]<-radialpos[i]+(labelsep-minsep)/2
    radialpos[i+1]<-radialpos[i+1]-(labelsep-minsep)/2
   }
  }
  xpos<-labelrad * radius * cos(radialpos[i])
  offset<-(radialpos[i] > pi && radialpos[i] < 2 * pi) * height
  ypos<-labelrad * radius * sin(radialpos[i]) * 2 * theta/pi +
   sin(radialpos[i]) * height
  text(xpos,ypos,labels[i],col=labelcol,
   adj=c(0.5,abs(0.5-sin(radialpos[i])/2)))
 }
 par(cex=oldcex,xpd=FALSE)
}

draw.tilted.sector<-function(x=0,y=0,edges=NA,radius=1,height=0.1,
 theta=pi/6,start=0,end=pi*2,border=par("fg"),col=par("bg"),
 explode=0,shade=0.8) {

 if(is.na(edges)) edges<-trunc(20*(end-start))
 angleinc<-(end-start)/edges
 angles<-c(seq(start,end,by=angleinc),end)
 viscurve<-(angles>=pi)&(angles<=2*pi)
 viscurvew<-(angles>=3*pi)
 nv<-length(angles)
 bisector<-(start+end)/2
 if(explode){
  # calculate the x and y offsets for the explode
  x<-x+cos(bisector)*explode
  y<-y+sin(bisector)*(1-sin(theta))*explode
 }
 if(shade>0 && shade<1){
  # calculate the shade color for the sides of the sector
  rgbcol<-col2rgb(col)
  shadecol<-rgb(shade*rgbcol[1]/255,shade*rgbcol[2]/255,
  shade*rgbcol[3]/255)
 }
 else shadecol<-col
 xp<-cos(angles) * radius + x
 # this is the top of the sector
 yp<-sin(angles) * 2 * theta/pi * radius + y
 if(start > 3*pi/2) {
  # the 'left' side will be visible in this quadrant
  if(explode > 0)
   # display the 'right' side just in case it goes beyond pi/2
   polygon(c(xp[nv],x,x,xp[nv],xp[nv]),c(yp[nv]-height,y-height,
    y+height,yp[nv]+height,yp[nv]+height),border=border,
    col=shadecol)
  # display the 'outside' of the sector
  polygon(c(xp[viscurve],rev(xp[viscurve])),c(yp[viscurve]-height,
   rev(yp[viscurve])+height),border=border,col=shadecol)
  # and any "wraparound" values
  polygon(c(xp[viscurvew],rev(xp[viscurvew])),c(yp[viscurvew]-height,
   rev(yp[viscurvew])+height),border=border,col=shadecol)
  if(explode > 0)
   # display the 'left' (front) side
   polygon(c(xp[1],x,x,xp[1],xp[1]),c(yp[1]-height,y-height,
    y+height,yp[1]+height,yp[1]),border=border,col=shadecol)
 }
 else {
  if(start > pi/2) {
   if(explode > 0) {
    polygon(c(xp[1],x,x,xp[1],xp[1]),c(yp[1]-height,y-height,
     y+height,yp[1]+height,yp[1]),border=border,
     col=shadecol)
    polygon(c(xp[nv],x,x,xp[nv],xp[nv]),c(yp[nv]-height,
     y-height,y+height,yp[nv]+height,yp[nv]+height),
     border=border,col=shadecol)
   }
   if(end > pi)
    polygon(c(xp[viscurve],rev(xp[viscurve])),c(yp[viscurve]-height,
     rev(yp[viscurve])+height),border=border,col=shadecol)
  }
  else {
   if(end > pi || start<2*pi)
    polygon(c(xp[viscurve],rev(xp[viscurve])),c(yp[viscurve]-height,
     rev(yp[viscurve])+height),border=border,
     col=shadecol)
   if(end > pi/2 && end < 3*pi/2 && explode > 0){
    polygon(c(xp[nv],x,x,xp[nv],xp[nv]),c(yp[nv]-height,
     y-height,y+height,yp[nv]+height,yp[nv]+height),
     border=border,col=shadecol)
   }
   if(explode > 0)
    polygon(c(xp[1],x,x,xp[1],xp[1]),c(yp[1]-height,y-height,
     y+height,yp[1]+height,yp[1]+height),border=border,
     col=shadecol)
  }
 }
 #display the top
 polygon(c(xp,x),c(yp+height,y+height),border=border,col=col)
 return(bisector)
}

pie3D<-function(x,edges=NA,radius=1,height=0.1,theta=pi/6, 
 start=0,border=par("fg"),col=NULL,labels=NULL,labelpos=NULL,
 labelcol=par("fg"),labelcex=1.5,labelrad=1.25,sector.order=NULL,
 explode=0,shade=0.8,mar=c(4,4,4,4),pty="s",...) {

 if(any(is.na(x))) x<-x[!is.na(x)]
 if(is.null(labels)) labels<-names(x)
 if(!is.numeric(x) || any(x < 0)) 
  stop("pie3D: x values must be positive numbers")
 # drop NAs
 oldmar<-par("mar")
 par(pty=pty,mar=mar,xpd=TRUE)
 x<-c(cumsum(c(0,x))/sum(x))*2*pi+start
 labelrad<-labelrad+explode
 lenx<-length(x)
 # add the last edge
 x[lenx]<-2*pi+start
 nsectors<-length(x)-1
 if(is.null(col)) col <- rainbow(nsectors)
 else if(length(col) < nsectors) col<-rep(col,length.out=nsectors)
 if(is.null(sector.order))
  # get the order of drawing sectors
  sector.order<-
   order(sin((x[2:lenx]+x[1:(lenx-1)])/2),decreasing=TRUE)
 bc<-rep(0,nsectors)
 # set up an empty plot
 plot(0,xlab="",ylab="",xlim=c(-1,1),ylim=c(-1,1),type="n",
  axes=FALSE,...)
 for(i in sector.order) {
  if(x[i] != x[i+1])
   bc[i]<-draw.tilted.sector(radius=radius,height=height, 
    theta=theta,start=x[i],end=x[i+1],edges=edges, 
    border=border,col=col[i],explode=explode,shade=shade)
  else bc[i]<-x[i]
 }
 if(!is.null(labels)) {
  if(!is.null(labelpos)) bc<-labelpos
  pie3D.labels(bc,radius=radius,height=height,theta=theta, 
   labels=labels,labelcol=labelcol,labelcex=labelcex,
   labelrad=labelrad)
 }
 par(mar=oldmar,xpd=FALSE,pty="m")
 invisible(bc)
}
