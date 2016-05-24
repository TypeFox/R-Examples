subnetworks <-
function(theta)
{
K = length(theta)
subnets = list()
for(k in 1:K)
{
	subnets[[k]] = list()
	m = theta[[k]]
		
	p=dim(m)[1]
	used = c()
	tempblock = c()  # vector of nodes in the current block
	indices = 1:(p+1)

	i=1
	snumber = 1
	iters = 1
	while((i<p+1)&(iters<p+1))
	{
		# find everything connected to i:
		tempblock = (1:p)==i
		tempblocksum = 0

		# loop to find all neighbors of neighbors of...
		while(tempblocksum!=sum(tempblock))
		{
			# save the size of the block at the beginning of this iteration
			tempblocksum = sum(tempblock) 
			# matrix of all neighbors of tempblock
			neighbors = m[tempblock,]!=0
			# update tempblock to include all new neighbors
			tempblock = colSums(rbind(neighbors,rep(0,p)))>0
		}
	
		# save the names of the features in tempblock:
		tb = names(tempblock)[tempblock]
		if(length(tb)>1)
		{
			subnets[[k]][[snumber]] = tb
			snumber = snumber+1
		}

		# mark the nodes as used:
		used = c(used,(1:p)[tempblock])
	
		# update i to lowest unused index:
		i = min(setdiff(indices,used))
		iters=iters+1
	}	#end loop
}
return(subnets)
}

