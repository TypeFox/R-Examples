jaccardmat <- function(res1,res2){
	if(res1@Number > 0 & res2@Number > 0){
		mat <- matrix(nrow=res1@Number,ncol=res2@Number)
		rownames(mat) <- paste("BC",1:res1@Number,sep="")
		colnames(mat) <- paste("BC",1:res2@Number,sep="")
		for(i in 1:res1@Number){
			for(j in 1:res2@Number){
				A <- which(res1@RowxNumber[,i]%*%t(res1@NumberxCol[i,])>0)
				B <- which(res2@RowxNumber[,j]%*%t(res2@NumberxCol[j,])>0)
				C <- intersect(A,B)
				mat[i,j] <- length(C)/(length(A)+length(B)-length(C))
			}
		}
	}
	else{
	mat <- matrix(0,ncol=1,nrow=1)
	}
	return(mat)
}
