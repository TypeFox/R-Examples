#   
#   Reduced model design matrix.
#
#   This function is called by the routine droptest.
#   Obtains the reduced model design matrix
#   for full model design matrix xmat and 
#   hypothesis matrix amat.  Checks to see
#   if amat has full row rank.
#   Uses Theorem 3.7.2 of H&M (1998, p. 186)
#
redmod<-function(xmat,amat){
  xmat=as.matrix(xmat)
  amat=rbind(amat)
  q<-length(amat[,1])
  p<-length(xmat[1,])
  temp<-qr(t(amat))
  if(temp$rank != q) {
    stop("redmod:  The hypothesis matrix is not full row rank.")
  } else {
    zed<-qr.qty(temp,t(xmat))
    redmod<-rbind(zed[(q+1):p,])
  }
  t(redmod)
}
