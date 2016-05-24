std.error<-function(x,na.rm) {
 vn<-function(x) return(sum(!is.na(x)))
 dimx<-dim(x)
 if(is.null(dimx)) {
  stderr<-sd(x,na.rm=TRUE)
  vnx<-vn(x)
 }
 else {
  if(is.data.frame(x)) {
   vnx<-unlist(sapply(x,vn))
   stderr<-unlist(sapply(x,sd,na.rm=TRUE))
  }
  else {
   vnx<-unlist(apply(x,2,vn))
   stderr<-unlist(apply(x,2,sd,na.rm=TRUE))
  }
 }
 return(stderr/sqrt(vnx))  
}
