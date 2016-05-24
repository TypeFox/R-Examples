"phi.range" <-
function(mat){
  cmat<-cor(mat)+diag(NA,ncol(mat))
  ma<-max(cmat,na.rm=TRUE)
  mi<-min(cmat,na.rm=TRUE)
  RET <- ma-mi
  RET
}

