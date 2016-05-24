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



#Aufbau der Hesse mit C
build.hess<-function(cm,im,deltas){
  n<-dim(cm)[1]
  nr3<-3*n
  dx=deltas$dx
  dy=deltas$dy
  dz=deltas$dz
  ds=deltas$ds
  hessian.mat<-matrix(data=0,ncol=nr3,nrow=nr3)
  out<-.C("buildHess",n=as.integer(n),contact.mat=as.integer(cm),hessian.mat=as.double(hessian.mat),interaction.mat=as.double(im),dx=as.double(dx),dy=as.double(dy),dz=as.double(dz),ds=as.double(ds),PACKAGE="BioPhysConnectoR")
  hessian.mat<-matrix(data=out$hessian.mat,nrow=nr3,ncol=nr3)
  return(hessian.mat)
 }

