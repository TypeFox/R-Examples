pdsXty <-
function(Xm,yv){
  xeig <- eigen(Xm,symmetric=TRUE)
  nze <- sum(xeig$val>xeig$val[1]*.Machine$double.eps)
  cpv <- (t(xeig$vec[,1:nze])*sqrt(xeig$val[1:nze]))%*%yv
}
