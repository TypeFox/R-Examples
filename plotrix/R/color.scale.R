color.scale<-function(x,cs1=c(0,1),cs2=c(0,1),cs3=c(0,1),alpha=1,
 extremes=NA,na.color=NA,xrange=NULL,color.spec="rgb") {

 if(diff(range(x,na.rm=TRUE))==0) x<-x/max(x,na.rm=TRUE)
 naxs<-is.na(x)
 if(!is.na(extremes[1])){
  # calculate the color ranges from the extremes - only for rgb
  colmat<-col2rgb(extremes)
  cs1<-colmat[1,]/255
  cs2<-colmat[2,]/255
  cs3<-colmat[3,]/255
  color_spec<-"rgb"
 }
 maxcs1<-ifelse(color.spec=="hcl",360,1)
 maxcs2<-ifelse(color.spec=="hcl",100,1)
 maxcs3<-ifelse(color.spec=="hcl",100,1)
 ncolors<-length(x)
 if(is.null(xrange)) {
  xrange<-range(x,na.rm=TRUE)
  drop.extremes<-FALSE
 }
 else {
  if(xrange[1] > min(x,na.rm=TRUE) || xrange[2] < max(x,na.rm=TRUE))
   stop("An explicit range for x must include the range of x values.")
  x<-c(xrange,x)
  drop.extremes=TRUE
 }
 ncs1<-length(cs1)
 if(ncs1>1) {
  cs1s<-rep(cs1[ncs1],ncolors)
  xstart<-xrange[1]
  xinc<-diff(xrange)/(ncs1-1)
  for(seg in 1:(ncs1-1)){
   segindex<-which((x >= xstart) & (x <= (xstart+xinc)))
   cs1s[segindex]<-rescale(x[segindex],cs1[c(seg,seg+1)])
   xstart<-xstart+xinc
  }
  if(min(cs1s,na.rm=TRUE) < 0 || max(cs1s,na.rm=TRUE) > maxcs1)
   cs1s<-rescale(cs1s,c(0,maxcs1))
 }
 else cs1s<-rep(cs1,ncolors)
 ncs2<-length(cs2)
 if(ncs2>1) {
  cs2s<-rep(cs2[ncs2],ncolors)
  xstart<-xrange[1]
  xinc<-diff(xrange)/(ncs2-1)
  for(seg in 1:(ncs2-1)){
   segindex<-which((x >= xstart) & (x <= (xstart+xinc)))
   cs2s[segindex]<-rescale(x[segindex],cs2[c(seg,seg+1)])
   xstart<-xstart+xinc
  }
  if(min(cs2s,na.rm=TRUE) < 0 || max(cs2s,na.rm=TRUE) > maxcs2)
   cs2s<-rescale(cs2s,c(0,maxcs2))
 }
 else cs2s<-rep(cs2,ncolors)
 ncs3<-length(cs3)
 if(ncs3>1) {
  cs3s<-rep(cs3[ncs3],ncolors)
  xstart<-xrange[1]
  xinc<-diff(xrange)/(ncs3-1)
  for(seg in 1:(ncs3-1)){
   segindex<-which((x >= xstart) & (x <= (xstart+xinc)))
   cs3s[segindex]<-rescale(x[segindex],cs3[c(seg,seg+1)])
   xstart<-xstart+xinc
  }
  if(min(cs3s,na.rm=TRUE) < 0 || max(cs3s,na.rm=TRUE) > maxcs3)
   cs3s<-rescale(cs3s,c(0,maxcs3))
 }
 else cs3s<-rep(cs3,ncolors)
 if(drop.extremes) {
  cs1s<-cs1s[-(1:2)]
  cs2s<-cs2s[-(1:2)]
  cs3s<-cs3s[-(1:2)]
 }
 xdim<-dim(x)
 colors<-do.call(color.spec,list(cs1s,cs2s,cs3s,alpha=alpha))
 if(!is.null(xdim)) colors<-matrix(colors,nrow=xdim[1])
 if(length(naxs)) colors[naxs]<-na.color
 return(colors)
}
