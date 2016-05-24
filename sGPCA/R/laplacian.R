laplacian <- function(dims){
	
	
	
	
	if(length(dims) > 3){
		stop("dimens must be a vector with 1 2 or 3 elements")
	}
	if(length(dims) < 1){
		stop("dimens must be a vector with 1 2 or 3 elements")
	}
	
	
	Q = 0
	
	
	if(length(dims) == 1){
		return_vectors <- .Call("onedlap_R",as.integer(dims))
		
		i = return_vectors[[2]][which(return_vectors[[3]]!=0)]
		j = return_vectors[[1]][which(return_vectors[[3]]!=0)]
		p = return_vectors[[3]][which(return_vectors[[3]]!=0)]
		
		i = i+1
		j = j+1
		
		Q <- sparseMatrix(i,j,x=as.numeric(p))
		
		
		
	}
	
	
	
	if(length(dims) == 2){	
		
		
		return_vectors <- .Call("twodlap_R",as.integer(dims[1]),as.integer(dims[2]))
		
		i = return_vectors[[2]][which(return_vectors[[3]]!=0)]
		j = return_vectors[[1]][which(return_vectors[[3]]!=0)]
		p = return_vectors[[3]][which(return_vectors[[3]]!=0)]
		
		i = i+1
		j = j+1
		
		Q <- sparseMatrix(i,j,x=as.numeric(p))
		
	}		
	
	
	if(length(dims) == 3){	
		
		
		return_vectors <- .Call("threedlap_R",as.integer(dims[1]),as.integer(dims[2]),as.integer(dims[3]))
		
		i = return_vectors[[2]][which(return_vectors[[3]]!=0)]
		j = return_vectors[[1]][which(return_vectors[[3]]!=0)]
		p = return_vectors[[3]][which(return_vectors[[3]]!=0)]
		
		i = i+1
		j = j+1
		
		Q <- sparseMatrix(i,j,x=as.numeric(p))
		
	}		
	
	
		
	return(Q)
	
	
}
