distLaplacian <- function(matrix,window = 2){
	Q = 0
	
	if(as.integer(window) < 2){
		stop("windowSize must be greater than or equal to 2")
	}
		
	
	if(class(matrix)== "dgCMatrix"){
		
		ir = matrix@i
		jc = matrix@p
		pr = matrix@x
		
		return_vectors <- .Call("lap_sparse_R",pr,as.integer(ir),as.integer(jc),as.integer(window),as.integer(0),dim(matrix)[1])
		
		i = return_vectors[[2]][which(return_vectors[[3]]!=0)]
		j = return_vectors[[1]][which(return_vectors[[3]]!=0)]
		p = return_vectors[[3]][which(return_vectors[[3]]!=0)]
		
		i = i+1
		j = j+1
		
		Q <- sparseMatrix(i,j,x=as.numeric(p))
	}else{
		
		
		return_vectors <- .Call("distance_lap_R",matrix,as.integer(window),as.integer(0),as.integer(dim(matrix)[1]))
		
		i = return_vectors[[2]][which(return_vectors[[3]]!=0)]
		j = return_vectors[[1]][which(return_vectors[[3]]!=0)]
		p = return_vectors[[3]][which(return_vectors[[3]]!=0)]
		
		i = i+1
		j = j+1
		
		Q <- sparseMatrix(i,j,x=as.numeric(p))
	}
	return(Q)
	
	
}
