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


 
get.svd<-function(hessian.mat, linpack=TRUE){
  s<-svd(hessian.mat,LINPACK=linpack)
  z<-matrix(data=0,nrow=length(s$d),ncol=2)
  z[,1]<-1:length(s$d)
  z[,2]<-s$d
 
 #Sortierung der Eigenwerte bei Beibehaltung der Indizes, z entspricht w
  z<-mat.sort(z,2)
  v<-s$v
  u<-s$u
  indx<-z[,1]
  ev<-z[,2]
  
  ret<-list(v=v,indx=indx,ev=ev,u=u)
  return(ret)
} 
