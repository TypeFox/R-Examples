# thigmophobe returns the direction (as 1|2|3|4 - see pos= in the text function) 
# _away_ from the nearest point where x and y are vectors of 2D coordinates

thigmophobe<-function(x,y) {
 # if x has at least two columns, split it
 if(missing(y)) {
  if(is.list(x) && length(x) >= 2) {
   y<-x[[2]]
   x<-x[[1]]
  }
  else {
   if(is.matrix(x) && dim(x)[2] >= 2) {
    y<-x[,2]
    x<-x[,1]
   }
   else
    stop("if y is missing, x must be a list with at least 2 columns")
  }
 }
 if(is.array(x)) x<-as.numeric(x)
 if(is.array(y)) y<-as.numeric(y)
 # get the current upper and lower limits of the plot
 plot.span<-par("usr")
 x.span<-plot.span[2] - plot.span[1]
 y.span<-plot.span[4] - plot.span[3]
 # if either axis is logarithmic, transform the values into logarithms
 if(par("xlog")) x<-log(x)
 if(par("ylog")) y<-log(y)
 # scale the values to the plot span
 # this avoids the numerically larger
 # axis dominating the distance measure
 x<-x/x.span
 y<-y/y.span
 # trash any names that may be attached to x or y
 names(x)<-names(y)<-NULL
 # get the distance matrix as a full matrix
 xy.dist<-as.matrix(dist(cbind(x,y)))
 lenx<-length(x)
 nearest.index<-rep(0,lenx)
 for(index in 1:lenx)
  nearest.index[index]<-as.numeric(names(which.min(xy.dist[-index,index])))
 # get the x and y differences for each point to the nearest point
 xdiff<-x - x[nearest.index]
 ydiff<-y - y[nearest.index]
 # first set the east/west direction
 dir.ew<-ifelse(xdiff > 0,4,2)
 # now do the north/south
 dir.ns<-ifelse(ydiff > 0,3,1)
 dir.away<-ifelse(abs(xdiff)>abs(ydiff),dir.ew,dir.ns)
 # set any congruent points to N/S labels or they'll overprint
 for(i in 1:lenx) {
  if(!xdiff[i] & !ydiff[i])
   dir.away[c(i,nearest.index[i])]<-c(1,3)
 }
 return(dir.away)
}

# thigmophobe.labels positions labels at points so that they
# are most distant from the nearest other point, where the
# points are described as x and y coordinates.

thigmophobe.labels<-function(x,y,labels=NULL,text.pos=NULL,...) {
 if(missing(x))
  stop("Usage: thigmophobe.labels(x,y,labels=1:length(x))")
 lenx<-length(x)
 # if x has at least two columns, split it
 if(missing(y)) {
  if(is.list(x) && lenx >= 2) {
   y<-x[[2]]
   x<-x[[1]]
  }
  else
   stop("if y is missing, x must be a list with at least 2 elements")
 }
 # check for NA or NaN
 validxy<-!(is.na(x) | is.na(y))
 if(is.null(labels)) labels<-1:lenx
 if(is.null(text.pos)) {
  if(lenx > 1) text.pos<-thigmophobe(x[validxy],y[validxy])
  else text.pos<-3
 }
 # allow labels to extend beyond the plot area
 par(xpd=TRUE)
 text(x[validxy],y[validxy],labels[validxy],pos=text.pos,...)
 # restore the clipping
 par(xpd=FALSE)
 invisible(list(x=x[validxy],y=y[validxy]))
}
