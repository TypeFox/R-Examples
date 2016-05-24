cluster.overplot<-function(x,y,away=NULL,tol=NULL,...) {
 if(missing(x)) stop("Usage: cluster.overplot(x,y,away=NULL,tol=NULL)")
 dimx<-dim(x)
 if(missing(y) && !is.null(dimx)) {
  y<-x[,2]
  x<-x[,1]
 }
 xlen<-length(x)
 if(xlen != length(y)) stop("x and y must be the same length.")
 if(dev.cur() == 1) plot(x,y,main="",xlab="",ylab="",axes=FALSE,type="n")
 if(is.null(away)) away<-c(strwidth("o"),5*strheight("o")/8)
 if(is.null(tol)) tol<-c(strwidth("o")/2,strheight("o")/2)
 away.x<-c(0,-away[1],away[1],0,0,-away[1],away[1],-away[1],away[1])
 away.y<-c(0,0,0,-away[2],away[2],-away[2],away[2],away[2],-away[2])
 flags<-1:xlen
 for(i in 1:xlen) {
  if(!is.na(flags[i])) {
   overplots<-abs(x-x[i]) <= tol[1] & abs(y - y[i]) <= tol[2]
   if(sum(overplots) > 1) {
    away.index<-1
    for(j in 1:xlen) {
     if(overplots[j]) {
      x[j]<-x[j]+away.x[away.index]
      y[j]<-y[j]+away.y[away.index]
      away.index<-away.index+1
      if(away.index > 9) away.index<-1
     }
    }
   }
  }
  flags[overplots] <- NA
 }
 return(list(x=x,y=y,...))
}
