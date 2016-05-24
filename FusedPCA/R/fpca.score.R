fpca.score <- function(A, maxsteps = 100, tol = 1e-3, K = 2, ridge = T, approx = F){
	
	# we need the diagonal elements equal to 0 
	diag(A) = 0
	
	# get the iso and non-iso lists
	iso.A = isolate(A)
	iso.seq = iso.A$isolate
	noniso.seq = iso.A$nonisolate
	
	# work on the non-iso part from now on
	A.noniso = A[noniso.seq, noniso.seq]
	
	# to get the laplacian matrix
	A.noniso.svd = svd(A.noniso)
	d.abs.order = sort(abs(A.noniso.svd$d), decreasing = T)
	leading.ind = apply(as.matrix(d.abs.order[1:K]), 1, function(x){which(abs(A.noniso.svd$d) == x)})
	#leading.ind = which(abs(A.noniso.svd$d) == d.abs.order[1:K])
	
	n.ind = length(leading.ind)
	
	n = length(noniso.seq)
	Ts = fused.trans(A.noniso)
	
	fused.whole = list()
	for(i in 1 : n.ind){
		
		# response vector in the equivalent regression case
	    reg.y = A.noniso.svd$u[, leading.ind[i]] * A.noniso.svd$d[leading.ind[i]]
	    
	    temp = fusedlasso.mod(y = reg.y, X = A.noniso, D = Ts, maxsteps = maxsteps, tol = tol, ridge = ridge, approx = approx)
	    fused.whole[[i]] = temp$beta
	}
	
	temp1 = rep(0, n.ind)
	for(i in 1: n.ind){
		temp1[i] = dim(fused.whole[[i]])[2]
	}
	ind.min = min(temp1)
	for(i in 1: n.ind){
		fused.whole[[i]] = fused.whole[[i]][, 1:ind.min]
	}
	final.whole = array(0, dim = c(n, n.ind, ind.min))
	final.matrix = array(0, dim = c(n, n.ind, ind.min))
	for(i in 1:ind.min){
		for(j in 1: n.ind){
			final.whole[, j, i] = fused.whole[[j]][,i]
		}
		final.matrix[, , i] = scale(final.whole[, , i], center = F)
	}
	
	final.r.matrix = array(0, dim = c(n, n.ind - 1, ind.min))
	for(i in 1:ind.min){
		for(j in 2: (n.ind - 1)){
			final.r.matrix[, (j-1), i] = (final.whole[, j, i])/(final.whole[, 1, i])
		}
		final.r.matrix[, , i] = scale(final.r.matrix[, , i], center = F)
	}
	final.whole = final.r.matrix
	class(final.whole) = 'FPCA-RoE'
	#temp = list(final.matrix = final.matrix, final.r.matrix = final.r.matrix)
	#class(temp) = 'FPCA-RoE'
	return(list(final.whole = final.whole, iso.seq = iso.seq, final.matrix = final.matrix))
}

