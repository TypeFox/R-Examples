fused.trans <- function(A){
	
	n = dim(A)[1]
	n.non.zero = sum(A[upper.tri(A)] != 0)
	T1 = matrix(0, nrow = n.non.zero, ncol = n)
	temp = 1
	for(i in 1: (n - 1)){
		for(j in (i + 1): n){
			if(A[i, j] == 1){
			    T1[temp, i] = 1
			    T1[temp, j] = -1
			    temp = temp + 1
			    
			} 
		}
	}
	
    return(T1)
}
