net.neighbors <-
function(theta,index)
{
	#index = (row.names(theta[[1]])==name)
	K = length(theta)
	p = dim(theta[[1]])[1]
	neighbors = list()
	for(k in 1:K)
	{
		neighbor.indices = (theta[[k]][index,]!=0)
		neighbor.indices[index] = FALSE
		neighbors[[k]] = row.names(theta[[1]])[neighbor.indices]
	}
	return(neighbors)
}

