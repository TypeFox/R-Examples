blkdiag <- function(A,B){
	
	dimA = dim(A)
	dimB = dim(B)
	C = array(0,dimA+dimB)
	C[1:dimA[1],1:dimA[2]] = A	
	C[dimA[1]+(1:dimB[1]),dimA[2]+(1:dimB[2])] = B
	return(C)
		
}