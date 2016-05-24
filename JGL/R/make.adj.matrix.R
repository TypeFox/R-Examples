make.adj.matrix <-
function(theta, separate=FALSE)
{
	K = length(theta)
	adj = list()
	if(separate)
	{
		for(k in 1:K)
		{
			adj[[k]] = (abs(theta[[k]])>1e-5)*1
		}
	}
	if(!separate)
	{
		adj = 0*theta[[1]]
		for(k in 1:K)
		{
			adj = adj+(abs(theta[[k]])>1e-5)*2^(k-1)
		}
	}
return(adj)
}

