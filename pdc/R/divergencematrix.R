

divergencematrix <-
function(codebooks, divergence)
{
	warning("Calls to divergencematrix(...) is deprecated!")
	return(divergenceMatrix(codebooks, divergence));
}

divergenceMatrix <-
function(codebooks, divergence)
{
	l <- dim(codebooks)[1]
	mt <- matrix(rep(0,l*l),l,)
	for (i in 1:l)
	{
		for (j in 1:l)
		{
			mt[i,j] <- divergence(codebooks[i,], codebooks[j,])
			mt[j,i] <- mt[i,j]	
		}
	}
	return(mt);
}
