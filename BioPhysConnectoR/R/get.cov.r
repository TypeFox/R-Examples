#############################################
#   This code is subject to the license as stated in DESCRIPTION.
#   Using this software implies accceptance of the license terms:
#    - GPL 2
#
#   (C) by F. Hoffgaard, P. Weil, and K. Hamacher in 2009.
#
#  keul(AT)bio.tu-darmstadt.de
#
#
#  http://www.kay-hamacher.de
#############################################



get.cov<-function(cm,im,deltas){
	n<-dim(cm)[1]
	hessian.mat<-build.hess(cm,im,deltas)
	out<-get.svd(hessian.mat)
	covelation.mat<-build.invhess(out)
	return(covelation.mat)
	}
