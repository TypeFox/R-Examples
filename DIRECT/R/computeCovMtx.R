computeCovMtx <-
function (sdWICluster, sdTSampling, sdResidual, n.time, n.rep, DATA.FORMAT=0)
{
	if (n.rep==1)
	{
		Smatrix = matrix (1, nrow=n.time, ncol=n.time)
		Imatrix = matrix (0, nrow=n.time, ncol=n.time)
		diag (Imatrix) = 1

		cov.matrix = sdWICluster^2 * Smatrix + sdTSampling^2 * Imatrix
	}
	else if (DATA.FORMAT==0)
	{
		# each entry is at least sdWICluster^2
		mtx.dim = n.time * n.rep
		cov.matrix = matrix (sdWICluster^2, nrow=mtx.dim, ncol=mtx.dim)
		
		# compute diagonal blocks
		diag.block = matrix (sdTSampling^2, nrow=n.rep, ncol=n.rep)
		diag (diag.block) = sdResidual^2 + sdTSampling^2
		
		for (t in 1:n.time)
			cov.matrix[n.rep*(t-1)+1:n.rep, n.rep*(t-1)+1:n.rep] = cov.matrix[n.rep*(t-1)+1:n.rep, n.rep*(t-1)+1:n.rep] + diag.block
	}
	else
	{
		Smatrix = matrix (1, nrow=n.time, ncol=n.time)
		Imatrix = matrix (0, nrow=n.time, ncol=n.time)
		diag (Imatrix) = 1
		diagonal.block = sdWICluster^2 * Smatrix + (sdTSampling^2 + sdResidual^2) * Imatrix
		offdiagonal.block = sdWICluster^2 * Smatrix + sdTSampling^2 * Imatrix
		
		cov.matrix = matrix (0, nrow=n.time*n.rep, ncol=n.time*n.rep)
		for (r in 1:n.rep)
		{
			cov.matrix[(n.time * (r-1) + 1:n.time), (n.time * (r-1) + 1:n.time)] = diagonal.block 
		}
		for (r1 in 1:(n.rep-1))
		{
			for (r2 in (r1+1):n.rep)
			{
				cov.matrix[(n.time * (r1-1) + 1:n.time), (n.time * (r2-1) + 1:n.time)] = offdiagonal.block 
				cov.matrix[(n.time * (r2-1) + 1:n.time), (n.time * (r1-1) + 1:n.time)] = offdiagonal.block 
			}
		}
	}
	
	return (cov.matrix)
	
}

