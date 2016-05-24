count.overplot<-function(x,y,tol=NULL,col=par("fg"),pch="1",...) {
 if(missing(x)) stop("Usage: count.overplot(x,y,tol=NULL,...)")
 dimx<-dim(x)
 if(missing(y)) {
  if(is.list(x) && names(x)[1] == "x") {
   y<-x[[2]]
   x<-x[[1]]
  }
  else {
   if(!is.null(dimx)) {
    y<-x[,2]
    x<-x[,1]
   }
  }
 }
 if(any(is.na(x) | is.na(y))) {
  indices<-!is.na(x) & !is.na(y)
  x<-x[indices]
  y<-y[indices]
 }
 xlim<-range(x)
 ylim<-range(y)
 xlen<-length(x)
 if(xlen != length(y)) stop("x and y must be the same length.")
 if(is.null(tol)) {
  if(dev.cur() == 1 ) plot(x,y,type="n",axes=FALSE,xlab="",ylab="")
  tol<-c(strwidth("o")/2,strheight("o")/2)
 }
 else {
  if (length(tol) == 1) tol <- rep(tol,2)
 }
 if(length(col) < xlen) col<-rep(col,xlen)
 if(length(pch) < xlen) pch<-rep(pch,xlen)
 flags<-1:xlen
 xsep<-ysep<-xdup<-ydup<-xydup<-sepcol<-seppch<-rep(0,xlen)
 nsep<-ndup<-0
 for(i in 1:xlen) {
  if(!is.na(flags[i])) {
   dups<-abs(x - x[i]) <= tol[1] & abs(y - y[i]) <= tol[2]
   ndups<-sum(dups)
   if(ndups > 1) {
    ndup<-ndup + 1
    xydup[ndup]<-ndups
    xdup[ndup]<-x[i]
    ydup[ndup]<-y[i]
   }
   else {
    nsep<-nsep + 1
    xsep[nsep]<-x[i]
    ysep[nsep]<-y[i]
    sepcol[nsep]<-col[i]
    seppch[nsep]<-pch[i]
   }
  }
  flags[dups]<-NA
 }
 plot(xsep[1:nsep],ysep[1:nsep],xlim=xlim,ylim=ylim,
  col=sepcol[1:nsep],pch=seppch[1:nsep],...)
 text(xdup[1:ndup],ydup[1:ndup],xydup[1:ndup],...)
}
