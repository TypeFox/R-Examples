fpca.nonscore <- function(A, maxsteps = 100, tol = 1e-3, normalised = T, K = 2, ridge = T, approx = FALSE){
	
	# we need the diagonal elements equal to 0 
	diag(A) = 0
	
	# get the iso and non-iso lists
	iso.A = isolate(A)
	iso.seq = iso.A$isolate
	noniso.seq = iso.A$nonisolate
	
	# work on the non-iso part from now on
	A.noniso = A[noniso.seq, noniso.seq]
	
	# to get the laplacian matrix
	L = laplacian(A = A.noniso, normalised = normalised)
	L.svd = svd(L)
	
	n = length(noniso.seq)
	Ts = fused.trans(A.noniso)
	
	temp.dim = dim(Ts)
	Ts = as.numeric(Ts)
	Ts = matrix(Ts, byrow = F, nrow = temp.dim[1])
	
	#if(K == 2) index = n
	#if(K > 2) index = ((n - K + 1): n)
	index = ((n - K + 1): (n - 1))
	n.ind = length(index)
	
	fused.whole = list()
	for(i in 1 : n.ind){
		
		# response vector in the equivalent regression case
	    reg.y = L.svd$u[, index[i]] * L.svd$d[index[i]]
	    
	    temp = fusedlasso.mod(y = reg.y, X = L, D = Ts, maxsteps = maxsteps, tol = tol, ridge = ridge, approx = approx)
	    fused.whole[[i]] = temp$beta
	}
	
	if(K == 2) final.whole = fused.whole[[1]]
	
	if(K > 2){
		temp1 = rep(0, n.ind)
		for(i in 1: n.ind){
			temp1[i] = dim(fused.whole[[i]])[2]
		}
		ind.min = min(temp1)
		for(i in 1: n.ind){
			fused.whole[[i]] = fused.whole[[i]][, 1:ind.min]
		}
		final.whole = array(0, dim = c(n, n.ind, ind.min))
		for(i in 1:ind.min){
			for(j in 1: n.ind){
				final.whole[, j, i] = fused.whole[[j]][,i]
			}
			final.whole[, , i] = scale(final.whole[, , i], center = T)
		}
	}
	class(final.whole) = 'FPCA'
	return(list(final.whole = final.whole, iso.seq = iso.seq))
}

