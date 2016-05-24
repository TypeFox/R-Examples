emptyspace<-function(x,y=NULL) {
 xlim<-par("usr")[1:2]
 ylim<-par("usr")[3:4]
 if(is.null(y)) {
  if(is.list(x)) {
   # trust me, I'm a list
   y<-as.vector(x[[2]])
   x<-as.vector(x[[1]])
  }
  else stop("both x and y values must be supplied")
 }
 xyna<-which(is.na(x) | is.na(y))
 if(any(xyna)) {
  # get rid of any NAs
  x<-x[-xyna]
  y<-y[-xyna]
 }
 # check for negative values and offset them out
 xoffset<-yoffset<-0
 if(any(x < 0)) {
  xoffset<-min(x)
  x<-x-xoffset
  xlim<-xlim-xoffset
 }
 if(any(y < 0)) {
  yoffset<-min(y)
  y<-y-yoffset
  ylim<-ylim-yoffset
 }
 # here begins Ray Brownrigg's code
 ox<-order(x)
 x<-x[ox]
 y<-y[ox]
 x<-c(xlim[1],x,xlim[2])
 y<-c(ylim[1],y,ylim[2])
 omaxa<-0
 halfxspan<-diff(xlim)/2
 halfyspan<-diff(ylim)/2
 # added braces here for clarity
 for(i in 1:(length(x)-1)) {
  for(j in (i+1):length(x)) {
   x1<-x[i]
   x2<-x[j]
   XX<-x2 - x1
   # if(XX > halfxspan) break
   yy<-c(ylim[1],y[(i+1):(j-1)],ylim[2])
   oyy<-order(yy)
   yy<-yy[oyy]
   dyy<-diff(yy)
   # whichdyy<-(dyy <= halfyspan) & (dyy >= 0.5*XX) & (dyy <= 2*XX)
   # the next line fixes the problem with very large empty spaces
   whichdyy<-(dyy >= 0.5*XX) & (dyy <= 2*XX)
   wy1<-yy[whichdyy]
   if(length(wy1) > 0) {
    wy2<-yy[(1:length(dyy))[whichdyy]+1]
    k<-which.max(dyy[whichdyy])
    maxa<-(x2-x1)*(wy2[k]-wy1[k])
    if(maxa > omaxa) {
     omaxa<-maxa
     mx1<-x1
     mx2<-x2
     my1<-wy1[k]
     my2<-wy2[k]
    }
   }
  }
 }
 # if offsets were used, remove them from the final values
 if(xoffset < 0) {
  mx1<-mx1+xoffset
  mx2<-mx2+xoffset
 }
 if(yoffset < 0) {
  my1<-my1+yoffset
  my2<-my2+yoffset
 }
 # return the center of the rectangle
 return(list(x=(mx1+mx2)/2,y=(my1+my2)/2))
}
