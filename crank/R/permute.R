permute<-function(x) {
 lenx<-length(x)
 if(lenx > 2) {
  nrows<-gamma(lenx+1)
  perms<-matrix(NA,nrow=nrows,ncol=lenx)
  blocklen<-gamma(lenx)
  perms[,1]<-rep(x,each=blocklen)
  for(block in 1:lenx)
   perms[((block-1)*blocklen+1):(block*blocklen),2:lenx]<-permute(x[-block])
  return(perms)
 }
 else return(rbind(x,rev(x)))
}
