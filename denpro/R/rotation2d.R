rotation2d<-function(dendat,alpha){
  Rx<-matrix(0,2,2)
  Rx[1,]<-c(cos(alpha),-sin(alpha))
  Rx[2,]<-c(sin(alpha),cos(alpha)) 
  detdat<-Rx%*%t(dendat)
  detdat<-t(detdat)
  return(detdat)
}

