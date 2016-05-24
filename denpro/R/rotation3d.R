rotation3d<-function(dendat,alpha,beta,gamma){
  Rx<-matrix(0,3,3)
  Rx[1,]<-c(1,0,0)
  Rx[2,]<-c(0,cos(alpha),-sin(alpha))
  Rx[3,]<-c(0,sin(alpha),cos(alpha))
  Ry<-matrix(0,3,3)
  Ry[1,]<-c(cos(beta),0,sin(beta))
  Ry[2,]<-c(0,1,0)
  Ry[3,]<-c(-sin(beta),0,cos(beta))
  Rz<-matrix(0,3,3)
  Rz[1,]<-c(cos(gamma),-sin(gamma),0)
  Rz[2,]<-c(sin(gamma),cos(gamma),0)
  Rz[3,]<-c(0,0,1)
 
  detdat<-Rx%*%Ry%*%Rz%*%t(dendat)
  detdat<-t(detdat)
  return(detdat)
}






