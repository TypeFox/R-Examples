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



get.bfacs<-function(covmat){
	if(dim(covmat)[1]!=dim(covmat)[2]){
		print("This is no square matrix. Please check the result.")
		}
	a<-matrix(diag(covmat),byrow=TRUE,ncol=3)
	hb<-rowSums(a)
    return(hb)
}
