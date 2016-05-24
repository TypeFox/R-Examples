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



mat.read<-function(file.name,ij=FALSE,sym=FALSE){
	if(ij){
		x<-matrix(data=scan(file.name),ncol=3,byrow=TRUE)
		m<-min(min(x[,1]),min(x[,2]))
		if(m==0){
			x[,1]<-x[,1]+1
			x[,2]<-x[,2]+1
		}
		if(sym){
			mat <- matrix(data = 0, nrow = max(x[,1:2]), ncol = max(x[,1:2]))
		}else{
			mat<-matrix(data=0,nrow=max(x[,1]),ncol=max(x[,2]))
		}
		for(i in 1:dim(x)[1]){
			mat[x[i,1],x[i,2]]<-x[i,3]
			if(sym){
				mat[x[i,2],x[i,1]]<-x[i,3]
				}
			}
		}else{
			mat<-as.matrix(read.table(file.name))
			}
	return(mat)
	}
