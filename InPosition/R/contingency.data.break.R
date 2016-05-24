contingency.data.break <- function(DATA,boot=FALSE){#,constrained=FALSE){
	DATA<-as.matrix(DATA)
	n.obs <- sum(DATA)
	n.cols <- ncol(DATA)
	n.rows <- nrow(DATA)	
	col.data <- rep(1:n.cols%x%rep.int(1,n.rows),as.vector(DATA))
	row.data <- rep(rep.int(1,n.cols)%x%1:n.rows,as.vector(DATA))

	if(boot){
		#This is for much more complex designs. It can happen later.
#		if(constrained){
#			ret.data <- t(apply(DATA,1,function(internal.data){return(rmultinom(1,sum(internal.data),prob=internal.data))}))			
#		}else{
			boot.samp <- sample(length(col.data),replace=TRUE)
			ret.data <- t(diag(n.rows)[row.data[boot.samp],])%*%diag(n.cols)[col.data[boot.samp],]		
#		}
	}else{
		#real.data <- t(diag(n.rows)[row.data,])%*%diag(n.cols)[col.data,]		
		ret.data <- t(diag(n.rows)[row.data,])%*%diag(n.cols)[sample(col.data),]			
	}
	colnames(ret.data) <- colnames(DATA)
	rownames(ret.data) <- rownames(DATA)
	return(ret.data)
}