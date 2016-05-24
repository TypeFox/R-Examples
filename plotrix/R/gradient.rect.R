gradient.rect<-function(xleft,ybottom,xright,ytop,reds,greens,blues, 
 col=NULL,nslices=50,gradient="x",border=par("fg")) {

 if(is.null(col)) col<-color.gradient(reds, greens, blues, nslices)
 else nslices<-length(col)
 nrect<-max(unlist(lapply(list(xleft,ybottom,xright,ytop),length)))
 if(nrect > 1) {
  if(length(xleft) < nrect) xleft<-rep(xleft,length.out=nrect)
  if(length(ybottom) < nrect) ybottom<-rep(ybottom,length.out=nrect)
  if(length(xright) < nrect) xright<-rep(xright,length.out=nrect)
  if(length(ytop) < nrect) ytop<-rep(ytop,length.out=nrect)
  for(i in 1:nrect)
   gradient.rect(xleft[i],ybottom[i],xright[i],ytop[i],
    reds,greens,blues,col,nslices,gradient,border=border)
 }
 else {
  if (gradient == "x") {
   xinc <- (xright - xleft)/nslices
   xlefts <- seq(xleft, xright - xinc, length = nslices)
   xrights <- xlefts + xinc
   rect(xlefts,ybottom,xrights,ytop,col=col,lty=0)
   rect(xlefts[1],ybottom,xrights[nslices],ytop,border=border)
  }
  else {
   yinc <- (ytop - ybottom)/nslices
   ybottoms <- seq(ybottom, ytop - yinc, length = nslices)
   ytops <- ybottoms + yinc
   rect(xleft,ybottoms,xright,ytops,col=col,lty=0)
   rect(xleft,ybottoms[1],xright,ytops[nslices],border=border)
  }
 }
 invisible(col)
}
