### function to generate cockroach graph

gen.cr <- function(n1){
	
	n = 4 * n1
	A = matrix(0, nrow = n, ncol = n)
	
	for(i in 1: n1){
		
		A[i, i+1] = 1
		A[i+1, i] = 1
		A[2*n1 + i, 2*n1 + i + 1] = 1
		A[2*n1 + i + 1, 2*n1 + i] = 1
	}
	
	for(i in (1 + n1): (2*n1 - 1)){
		
		A[i, i+1] = 1
		A[i+1, i] = 1
		A[i, 2*n1 + i] = 1
		A[2*n1 + i, i] = 1
		A[2*n1 + i, 2*n1 + i + 1] = 1
		A[2*n1 + i + 1, 2*n1 + i] = 1
	}

	A[2*n1, 4*n1] = 1
	A[4*n1, 2*n1] = 1

	return(A)
	
}

