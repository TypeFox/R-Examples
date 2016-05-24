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



#Berechnung der Inversen Hesse mit C
 build.invhess<-function(svd_obj, singularity=6){
  ev<-svd_obj$ev
  nr3<-length(ev)
  precision<-10^(-9)
  ss<-length(which(ev<precision))
  if(!is.null(singularity)){
	if(ss!=singularity){
		cat("Found ",ss,"eigenvalues smaller than",precision,", and not",singularity, "use ",singularity,". Please check the result.\n")
	}
  }else{
	singularity<-ss
  }
  cov.mat<-matrix(data=0,ncol=nr3,nrow=nr3)
  out<-.C("invHess",nr3=as.integer(nr3),cov.mat=as.double(cov.mat),v=as.double(svd_obj$v),indx=as.integer(svd_obj$indx),ev=as.double(ev),singularity=as.integer(singularity),u=as.double(svd_obj$u),PACKAGE="BioPhysConnectoR")
  cov.mat<-matrix(data=out$cov.mat,nrow=nr3,ncol=nr3)
  return(cov.mat)
}
