hexagon<-function(x,y,unitcell=1,col=NA,border="black") {
 polygon(c(x,x,x+unitcell/2,x+unitcell,x+unitcell,x+unitcell/2),
  c(y+unitcell*0.125,y+unitcell*0.875,y+unitcell*1.125,y+unitcell*0.875,
    y+unitcell*0.125,y-unitcell*0.125),col=col,border=border)
}

fill.corner<-function(x,nrow,ncol,na.value=NA) {
 xlen<-length(x)
 ncells<-ifelse(nrow*ncol < xlen,nrow*ncol,xlen)
 newmat<-matrix(na.value,nrow=nrow,ncol=ncol)
 xside<-1
 while(xside*xside < ncells) xside<-xside+1
 row=1
 col=1
 for(xindex in 1:ncells) {
  newmat[row,col]<-x[xindex]
  if(row == xside) {
   col<-col+1
   row<-1
  }
  else row<-row+1
 }
 return(newmat)
}

color2D.matplot<-function(x,cs1=c(0,1),cs2=c(0,1),cs3=c(0,1),
 extremes=NA,cellcolors=NA,show.legend=FALSE,nslices=10,xlab="Column",
 ylab="Row",do.hex=FALSE,axes=TRUE,show.values=FALSE,vcol=NA,vcex=1,
 border="black",na.color=NA,xrange=NULL,color.spec="rgb",yrev=TRUE,
 xat=NULL,yat=NULL,Hinton=FALSE,...) {
 
 if(diff(range(x,na.rm=TRUE))==0) x<-x/max(x,na.rm=TRUE)
 if(is.matrix(x) || is.data.frame(x)) {
  xdim<-dim(x)
  if(is.data.frame(x)) x<-unlist(x)
  else x<-as.vector(x)
  oldpar<-par("xaxs","yaxs","xpd","mar")
  par(xaxs="i",yaxs="i")
  if(do.hex) par(mar=c(5,4,4,4))
  plot(c(0,xdim[2]),c(0,xdim[1]),xlab=xlab,ylab=ylab,type="n",axes=FALSE,...)
  oldpar$usr<-par("usr")
  if(!do.hex) {
   box()
   pos<-0
  }
  else pos<- -0.3
  if(axes) {
   if(is.null(xat)) xat<-pretty(0:xdim[2])[-1]
   axis(1,at=xat-0.5,labels=xat,pos=pos)
   if(is.null(yat)) yat<-pretty(0:xdim[1])[-1]
   axis(2,at=xdim[1]-yat+0.5,labels=yat)
  }
  if(all(is.na(cellcolors))) {
   if(Hinton) {
    if(is.na(extremes[1])) extremes<-c("black","white")
    cellcolors<-extremes[(x > 0) + 1]
   }
   else cellcolors<-color.scale(x,cs1,cs2,cs3,extremes=extremes,
    na.color=na.color,color.spec=color.spec)
  }
  # this sets the color for overprinted text to black or white
  # depending upon what color will be the background for the text
  if(is.na(vcol))
   vcol<-ifelse(colSums(col2rgb(cellcolors)*c(1,1.4,0.6))<350,"white","black")
  # if it's a Hinton diagram,cellsize = x, rescaling to 0,1 if necessary
  if(Hinton) {
   if(any(x < 0 | x > 1)) cellsize<-matrix(rescale(abs(x),c(0,1)),nrow=xdim[1])
  }
  else cellsize<-matrix(1,nrow=xdim[1],ncol=xdim[2])
  # start from the top left - isomorphic with the matrix layout
  if(do.hex) {
   par(xpd=TRUE)
   offset<-0
   if(length(border) < xdim[1]*xdim[2])
    border<-rep(border,length.out=xdim[1]*xdim[2])
   for(row in 1:xdim[1]) {
    for(column in 0:(xdim[2]-1)) {
     hexagon(column+offset,xdim[1]-row,unitcell=cellsize[row,column+1],
      col=cellcolors[row+xdim[1]*column],
      border=border[row+xdim[1]*column])
     if(show.values)
      text(column+offset+0.5,xdim[1]-row+0.5,x[row+column*xdim[1]],
       col=vcol[row+xdim[1]*column],cex=vcex)
    }
    offset<-ifelse(offset,0,0.5)
   }
   par(xpd=FALSE)
  }
  else {
   if(Hinton) inset<-(1-cellsize)/2
   else inset<-0
   if(yrev) {
    y0<-rep(seq(xdim[1]-1,0,by=-1),xdim[2])+inset
    y1<-rep(seq(xdim[1],1,by=-1),xdim[2])-inset
   }
   else {
    y0<-rep(0:(xdim[1]-1),xdim[2])+inset
    y1<-rep(1:xdim[1],xdim[2])-inset
   }
   rect(sort(rep((1:xdim[2])-1,xdim[1]))+inset,y0,
    sort(rep(1:xdim[2],xdim[1]))-inset,y1,
    col=cellcolors,border=border)
   if(show.values) {
    if(yrev) texty<-rep(seq(xdim[1]-0.5,0,by=-1),xdim[2])
    else texty<-rep(seq(0.5,xdim[1]-0.5,by=1),xdim[2])
    text(sort(rep((1:xdim[2])-0.5,xdim[1])),texty,
     round(x,show.values),col=vcol,cex=vcex)
   }
  }
  naxs<-which(is.na(x))
  xy<-par("usr")
  plot.din<-par("din")
  plot.pin<-par("pin")
  bottom.gap<-(xy[3]-xy[4])*(plot.din[2]-plot.pin[2])/(2*plot.pin[2])
  grx1<-xy[1]
  gry1<-bottom.gap*0.95
  grx2<-xy[1]+(xy[2]-xy[1])/4
  gry2<-bottom.gap*0.8
  if(length(cellcolors) > 1) {
   colmat<-col2rgb(c(cellcolors[which.min(x)],cellcolors[which.max(x)]))
   cs1<-colmat[1,]/255
   cs2<-colmat[2,]/255
   cs3<-colmat[3,]/255
   color.spec<-"rgb"
  }
  rect.col<-color.scale(1:nslices,cs1,cs2,cs3,color.spec=color.spec)
  if(show.legend)
   color.legend(grx1,gry1,grx2,gry2,round(range(x,na.rm=TRUE),show.legend),
    rect.col=rect.col)
  par(oldpar)
 }
 else cat("x must be a data frame or matrix\n")
}
