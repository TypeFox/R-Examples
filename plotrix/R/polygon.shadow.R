polygon.shadow<-function(x,y=NULL,offset=NA,inflate=NA,
 col=c("#ffffff","#cccccc")) {
 
 if(is.null(y)) {
  if(is.null(x$y)) stop("both x and y coordinates must be given")
  y<-x$y
  x<-x$x
 }
 xcenter<-mean(x)
 ycenter<-mean(y)
 if(is.na(offset[1])) offset<-c((max(x)-min(x))/20,-(max(y)-min(y))/20)
 if(length(col) == 2) col<-smoothColors(col[1],8,col[2])
 if(is.na(inflate[1])) {
  multx<-seq(offset[1],0,length=10)
  multy<-seq(offset[2],0,length=10)
 }
 else {
  multx<-seq(inflate[1],0,length=10)
  multy<-seq(inflate[2],0,length=10)
 }
 for(i in 1:10)
  polygon(x+multx[i]*(x-xcenter)+offset[1],y+multy[i]*(y-ycenter)+offset[2],
   col=col[i],border=NA)
}
