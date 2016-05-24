## function to generate degree-corrected block model
gen.dcbm <- function(n, theta.in, theta.bw, theta, K, seed){
	
	if(K == 2){
		set.seed(seed)
		
		# block for community 1 to 1
		A11 = matrix(0, ncol = n, nrow = n)
		for(i in 1:(n-1)){
			for(j in (i+1):n){
				A11[i, j] = rbinom(1, 1, theta.in*theta[i]*theta[j])
			}
		}
		
		# block for community 2 to 2
		A22 = matrix(0, ncol = n, nrow = n)
		for(i in 1:(n-1)){
			for(j in (i+1):n){
				A22[i, j] = rbinom(1, 1, theta.in*theta[n+i]*theta[n+j])
			}
		}
		
		# blcok for community 1 to 2
		A12 = matrix(0, ncol = n, nrow = n)
		for(i in 1:n){
			for(j in 1:n){
				A12[i, j] = rbinom(1, 1, theta.bw*theta[i]*theta[n+j])
			}
		}
		A1 = cbind(A11, A12)
		A2 = cbind(A12, A22)
		A = rbind(A1, A2)
		
		ind = lower.tri(A)
		A[ind] = t(A)[ind]
		diag(A) = rep(0, dim(A)[1])
	
	}
	
	if(K == 4){
		set.seed(seed)
		
		# block for community 1 to 1
		A11 = matrix(0, ncol = n, nrow = n)
		for(i in 1:(n-1)){
			for(j in (i+1):n){
				A11[i, j] = rbinom(1, 1, theta.in*theta[i]*theta[j])
			}
		}
		
		# block for community 2 to 2
		A22 = matrix(0, ncol = n, nrow = n)
		for(i in 1:(n-1)){
			for(j in (i+1):n){
				A22[i, j] = rbinom(1, 1, theta.in*theta[n+i]*theta[n+j])
			}
		}
		
		# block for community 3 to 3
		A33 = matrix(0, ncol = n, nrow = n)
		for(i in 1:(n-1)){
			for(j in (i+1):n){
				A33[i, j] = rbinom(1, 1, theta.in*theta[2*n+i]*theta[2*n+j])
			}
		}
		
		# block for community 4 to 4
		A44 = matrix(0, ncol = n, nrow = n)
		for(i in 1:(n-1)){
			for(j in (i+1):n){
				A44[i, j] = rbinom(1, 1, theta.in*theta[3*n+i]*theta[3*n+j])
			}
		}
		
		# blcok for community 1 to 2
		A12 = matrix(0, ncol = n, nrow = n)
		for(i in 1:n){
			for(j in 1:n){
				A12[i, j] = rbinom(1, 1, theta.bw*theta[i]*theta[n+j])
			}
		}
		
		# blcok for community 1 to 3
		A13 = matrix(0, ncol = n, nrow = n)
		for(i in 1:n){
			for(j in 1:n){
				A13[i, j] = rbinom(1, 1, theta.bw*theta[i]*theta[2*n+j])
			}
		}
		
		# blcok for community 1 to 4
		A14 = matrix(0, ncol = n, nrow = n)
		for(i in 1:n){
			for(j in 1:n){
				A14[i, j] = rbinom(1, 1, theta.bw*theta[i]*theta[3*n+j])
			}
		}
		
		# blcok for community 2 to 3
		A23 = matrix(0, ncol = n, nrow = n)
		for(i in 1:n){
			for(j in 1:n){
				A23[i, j] = rbinom(1, 1, theta.bw*theta[n+i]*theta[2*n+j])
			}
		}
		
		# blcok for community 2 to 4
		A24 = matrix(0, ncol = n, nrow = n)
		for(i in 1:n){
			for(j in 1:n){
				A24[i, j] = rbinom(1, 1, theta.bw*theta[n+i]*theta[3*n+j])
			}
		}
		
		# blcok for community 3 to 4
		A34 = matrix(0, ncol = n, nrow = n)
		for(i in 1:n){
			for(j in 1:n){
				A34[i, j] = rbinom(1, 1, theta.bw*theta[2*n+i]*theta[3*n+j])
			}
		}
		
		A1 = cbind(A11, A12, A13, A14)
		A2 = cbind(A12, A22, A23, A23)
		A3 = cbind(A13, A23, A33, A34)
		A4 = cbind(A14, A24, A34, A44)
		A = rbind(A1, A2, A3, A4)
		
		ind = lower.tri(A)
		A[ind] = t(A)[ind]
		diag(A) = rep(0, dim(A)[1])
	
	}

	
	
	
	return(A)
	
}