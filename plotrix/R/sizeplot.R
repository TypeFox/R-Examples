sizeplot<-function(x,y,scale=1,pow=0.5,powscale=TRUE,size=c(1,4),add=FALSE,...) {

 pair <- function(x, y) paste(x, y, sep = "/")
 unpair <- function(x) t(sapply(strsplit(x, "/"), as.numeric))
 f<-factor(pair(x, y))
 n<-table(f)
 if(min(n) == max(n)) {
  warning("all points repeated equally (why use sizeplot?)")
  if(powscale) psize <- rep(length(f), scale)
  else psize <- rep(length(f), size[1])
 }
 else {
  if(powscale) psize<-scale*n^pow
  else psize<-size[1]+(size[2]-size[1])*((n-min(n))/(max(n)-min(n)))
 }
 newpts<-unpair(levels(f))
 if(!add) plot(newpts[,1],newpts[,2],cex=psize,type="p",...)
 else points(newpts[,1],newpts[,2],cex=psize,...)
}
