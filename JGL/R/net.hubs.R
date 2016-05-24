net.hubs <-
function(theta,nhubs=10)
{
	degree = net.degree(theta)
	K = length(degree)
	hubs = list()
	for(k in 1:K)
	{
		order = order(degree[[k]],decreasing=TRUE)
		hubs[[k]] = degree[[k]][order[1:nhubs]]
	}
	return(hubs)
}

