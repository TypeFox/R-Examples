find_overplots<-function(x,y,tol) {
 xlen<-length(x)
 flags<-1:xlen
 xsep<-ysep<-xdup<-ydup<-xydup<-rep(NA,xlen)
 nsep<-ndup<-0
 for(i in 1:xlen) {
  if(!is.na(flags[i])) {
   dups<-abs(x - x[i]) <= tol[1] & abs(y - y[i]) <= tol[2]
   ndups<-sum(dups,na.rm=TRUE)
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
   }
  }
  flags[dups]<-NA
 }
 return(list(xsep=xsep,ysep=ysep,xdup=xdup,ydup=ydup,xydup=xydup))
}

display.overplots<-function(x,y,tol=NULL,how=c("count","cluster","size"),
 xlim=NULL,ylim=NULL,col=rep(par("fg"),2),pch=c("1",1),spc=NULL,...) {

 if(missing(x)) stop("display.overplots must have xy coordinates")
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
 if(is.null(xlim)) xlim<-range(x)
 if(is.null(ylim)) ylim<-range(y)
 xlen<-length(x)
 if(xlen != length(y)) stop("x and y must be the same length.")
 plot(x,y,type="n",axes=FALSE,xlab="",ylab="")
 xylim<-par("usr")
 if(is.null(tol)) tol<-c(diff(xylim[1:2]/100),diff(xylim[3:4]/100))
 else if(length(tol) == 1) tol <- rep(tol,2)
 if(length(col) < xlen) col<-rep(col,xlen)
 if(length(pch) < xlen) pch<-rep(pch,xlen)
 if(is.null(spc)) spc<-c(diff(xylim[1:2])/100,diff(xylim[3:4]/100))
 xy<-find_overplots(x,y,tol)
 plot(xy$xsep,xy$ysep,xlim=xlim,ylim=ylim,col=col[1],pch=pch[1],...)
 if(how[1] == "count") text(xy$xdup,xy$ydup,xy$xydup,col=col[2],...)
 if(how[1] == "cluster") {
  xshifts<-c(0,-spc[1],spc[1],0,0,-spc[1],spc[1],-spc[1],spc[1])
  yshifts<-c(0,0,0,-spc[2],spc[2],-spc[2],spc[2],spc[2],-spc[2])
  for(dup in 1:sum(!is.na(xy$xdup))) {
   for(ndups in 1:min(c(xy$xydup[dup],9))) {
    points(xy$xdup[dup]+xshifts[ndups],xy$ydup[dup]+yshifts[ndups],
     pch=pch[2],col=col[2],...)
   }
   if(xy$xydup[dup] > 9) points(xy$xdup[dup],xy$ydup[dup],pch=4,cex=2)
  }
 }
 if(how[1] == "size")
  points(xy$xdup,xy$ydup,pch=pch[2],col=col[2],cex=sqrt(xy$xydup))
}
