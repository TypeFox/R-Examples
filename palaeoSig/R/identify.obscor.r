identify.obscor<-function(x, labels, ...){
  if(missing(labels))labels<-rownames(x$ob$x)
  identify(x$ob$x[,1:2],labels=labels,...)
}
