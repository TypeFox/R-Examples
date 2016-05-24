net.degree <-
function(theta)
{
K = length(theta)
degree = list()
for(k in 1:K)
{
	degree[[k]] = (rowSums(abs(theta[[k]])>1e-5)-1)
}

return(degree)
}

