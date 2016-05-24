color.scale.lines<-function(x,y,reds,greens,blues,col=NA,colvar=NA,...) {

 if(is.list(x) && missing(y)) {
  y<-x$y
  x<-x$x
 }
 lenx<-length(x)
 leny<-length(y)
 nseg<-lenx-1
 if(lenx != leny) {
  warning("x and y are different lengths - some values will be ignored")
  if(lenx>leny) nseg<-leny-1
 }
 if(is.na(col[1])) {
  # if colors are not supplied
  if(is.na(colvar[1]))
   # if a separate variable is not supplied, use running average of y pairs
   lcol<-color.scale((y[1:nseg]+y[2:(nseg+1)])/2,reds,greens,blues)
  # assume that colvar refers to the intervals, not the x/y points
  else lcol<-color.scale(colvar,reds,greens,blues)
 }
 else lcol=col
 segments(x[1:nseg],y[1:nseg],x[2:(nseg+1)],y[2:(nseg+1)],col=lcol,...)
}
