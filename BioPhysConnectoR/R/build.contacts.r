#############################################
#   This code is subject to the license as stated in DESCRIPTION.
#   Using this software implies acceptance of the license terms:
#    - GPL 2
#
#   (C) by F. Hoffgaard, P. Weil, and K. Hamacher in 2009.
#
#  keul(AT)bio.tu-darmstadt.de
#
#
#  http://www.kay-hamacher.de
#############################################



#Aufbau der Distanzvektoren und der Kontaktmatrix in C
build.contacts<-function(n,cuts,xyz){
  contact.mat<-matrix(data=0,nrow=n,ncol=n)
  x<-xyz[,1]
  y<-xyz[,2]
  z<-xyz[,3]
  dx<-matrix(data=0,nrow=n,ncol=n)
  dy<-matrix(data=0,nrow=n,ncol=n)
  dz<-matrix(data=0,nrow=n,ncol=n)
  ds<-matrix(data=0,nrow=n,ncol=n)

out<-.C("contactDist",contact.mat=as.integer(contact.mat),dx=as.double(dx),dy=as.double(dy),dz=as.double(dz),ds=as.double(ds),n=as.integer(n),cuts=as.double(cuts),x=as.double(x),y=as.double(y),z=as.double(z),PACKAGE="BioPhysConnectoR")


  contact.mat<-matrix(data=out$contact.mat,nrow=n)
  dx<-matrix(data=out$dx,nrow=n)
  dy<-matrix(data=out$dy,nrow=n)
  dz<-matrix(data=out$dz,nrow=n)
  ds<-matrix(data=out$ds,nrow=n)
  count<-sum(contact.mat)/2
  delta<-list(dx=dx,dy=dy,dz=dz,ds=ds)
  ret<-list(cm=contact.mat,deltas=delta,cnr=count)

  return(ret)
}
