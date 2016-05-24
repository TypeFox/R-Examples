

divergencematrix.multichannel <-
function(codebooks, divergence, num.channels)
{
	warning("Calls to divergencematrix.multichannel(...) is deprecated!")
	return(divergenceMatrixMultichannel(codebooks, divergence, num.channels));
}

divergenceMatrixMultichannel <-
function(codebooks, divergence, num.channels)
{
	l <- dim(codebooks)[1]/num.channels
	mt <- matrix(rep(0,l*l),l,)
	for (i in 1:l)
	{
		for (j in 1:l)
		{
			mt[i,j] <- 0
			for (k in 1:num.channels) { 
				mt[i,j] <- mt[i,j] + (divergence(codebooks[(i-1)*num.channels+k,], codebooks[(j-1)*num.channels+k,]))**2
			}
			mt[i,j] <- sqrt(mt[i,j])
			mt[j,i] <- mt[i,j]	
		}
	}
	return(mt);
}
