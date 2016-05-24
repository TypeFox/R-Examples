statvc4levelt3 <-
function(y,cluster,m1,n1,m2,n2,m3,n3,p,weight){

# y = Matrix Nxp (N obs, p=dimension) of (maybe standardized; signs or ranks) data 
# cluster = List of indices (each element gives the indices of the obs in a level 3 cluster nested within a level 2 cluster nested within a level 1 cluster). From the "preparecluster" function.
# N = Total number of observations
# p = Dimension
# m1 = Number of level 1 clusters
# n1 = Number of observations in each level 1 cluster (vector m1X1)
# m2 = Number of level 2 clusters within each level 1 cluster (vector m1X1)
# n2 = Number of observations in each level 2 cluster within each level 1 cluster (list of size m1)
# m3 = Number of level 3 clusters nested within each combination of level 1 and level 2 clusters (list of size m1)
# n3 = Number of observations in each level 3 cluster nested within each combination of level 1 and level 2 clusters (list of size m1)


mat=matrix(0,p,p)
w=0

for(i in 1:m1)
	{
	for(j in 1:m2[i])
		{
		clusterij=cluster[[i]][[j]]
		m1ij=m3[[i]][j]
		n1ij=n3[[i]][[j]]		
		matij=statvc2level(y,clusterij,m1ij,n1ij,p,weight)
		mat=mat+matij[[1]]
		w=w+matij[[2]]
		}
	}

list(mat,w)
}
