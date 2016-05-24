tensorizza <-
function  (dif, mat) {
  dif <- as.matrix(dif)
  p <- ncol(dif)
  n <- nrow(dif)
  
  matrep <- array(mat, c(p,p,n))
  difrep<-array(dif,c(n,p,p))
  dift<-array(dif,c(n,p,p))
  matrep<-aperm(matrep,c(3,1,2))
  dift<-aperm(dift,c(1,3,2))
  res<-difrep*matrep*dift
  resfin<-apply(res,1,sum)
}

