gen.sbm <- function(n, theta.in, theta.bw, K, seed){
	
	if(K == 2){
		set.seed(seed)
		A1 = matrix(rbinom(n^2, 1, theta.in), nrow = n)
		set.seed(seed + 1)
		A2 = matrix(rbinom(n^2, 1, theta.bw), nrow = n)
        A = rbind(cbind(A1, A2), cbind(t(A2), A1))

		
	}
	
	if(K == 3){
		set.seed(seed)
		A1 = matrix(rbinom(n^2, 1, theta.in), nrow = n)
		set.seed(seed + 1)
		A5 = matrix(rbinom(n^2, 1, theta.in), nrow = n)
		set.seed(seed + 2)
		A6 = matrix(rbinom(n^2, 1, theta.in), nrow = n)
		set.seed(seed + 3)
		A2 = matrix(rbinom(n^2, 1, theta.bw), nrow = n)
		set.seed(seed + 4)
		A3 = matrix(rbinom(n^2, 1, theta.bw), nrow = n)
		set.seed(seed + 5)
		A4 = matrix(rbinom(n^2, 1, theta.bw), nrow = n)
		A = rbind(cbind(A1, A2, A3), cbind(t(A2), A5, A4), cbind(t(A3), t(A4), A6))

	}
	
	if(K == 4){
		set.seed(seed)
		A1 = matrix(rbinom(n^2, 1, theta.in), nrow = n)
		set.seed(seed + 1)
		A5 = matrix(rbinom(n^2, 1, theta.in), nrow = n)
		set.seed(seed + 2)
		A6 = matrix(rbinom(n^2, 1, theta.in), nrow = n)
		set.seed(seed + 3)
		A7 = matrix(rbinom(n^2, 1, theta.in), nrow = n)
		set.seed(seed + 7)
		A2 = matrix(rbinom(n^2, 1, theta.bw), nrow = n)
		set.seed(seed + 4)
		A3 = matrix(rbinom(n^2, 1, theta.bw), nrow = n)
		set.seed(seed + 5)
		A4 = matrix(rbinom(n^2, 1, theta.bw), nrow = n)
		set.seed(seed + 6)
		A8 = matrix(rbinom(n^2, 1, theta.bw), nrow = n)
		set.seed(seed + 8)
		A9 = matrix(rbinom(n^2, 1, theta.bw), nrow = n)
		set.seed(seed + 9)
		A10 = matrix(rbinom(n^2, 1, theta.bw), nrow = n)
		set.seed(seed + 10)
		A = rbind(cbind(A1, A2, A3, A8), cbind(t(A2), A5, A4, A9), cbind(t(A3), t(A4), A6, A10), cbind(t(A8), t(A9), t(A10), A7))

	}
	
	ind = lower.tri(A)
	A[ind] = t(A)[ind]
	diag(A) = rep(0, dim(A)[1])
	
	return(A)
	
}
