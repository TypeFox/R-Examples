mat.write<-function(mat,file.name,ij=FALSE,sym=FALSE,sparse=FALSE,formatted=TRUE){
	if(!ij){
		write.table(mat,file=file.name,quote=FALSE,row.names=FALSE,col.names=FALSE)
	}else{
		d<-dim(mat)
		# for output of non-zero values only
		if(sparse){
			ind<-which(mat!=0,arr.ind=TRUE)
			mat2<-cbind(ind,mat[ind])
		}else{
			j<-rep(1:d[1],d[2])
			i<-sort(j)
			ind<-cbind(i,j)
			mat2<-cbind(ind,mat[ind])
		}
		# if matrix is symmetric and only the upper triangle and the diagonal is needed
		if(sym){
			mat2<-mat2[which(mat[,1]<=mat[,2]),]
		}
		# if the output should be formatted (e.g. for gnuplot use)
		if(formatted){
			ind<-mat2[1,1]
			for(i in 1:dim(mat2)[1]){
				if(ind!=mat2[i,1]){
					cat("\n",file=file.name,append=TRUE)
					ind<-mat2[i,1]
				}
				cat(mat2[i,],"\n",sep=" ",file=file.name,append=TRUE)
			}
		}else{
			write.table(mat2,file=file.name,quote=FALSE,row.names=FALSE,col.names=FALSE)
		}
	}
}
