testvc2level <-
function(y,cluster,npermut,N,p,smat,m1,n1,weight){

# y = Matrix Nxp (N obs, p=dimension) of (maybe standardized; signs or ranks) data 
# cluster = List of indices (each element gives the indices of the obs in a level 1 cluster). From the "preparecluster" function.
# npermut = Number of permutations for the test
# N = Total number of observations
# p = Dimension
# smat = Transformation matrix to achieve affine-equivariance (from the "ystand" function)
# m1 = Number of level 1 clusters
# n1 =  Number of observations in each level 1 cluster (vector m1X1)
# weight = weights used to compute the dependence matrix
#		"observation" gives equal weight to each observation
#		"pair" gives equal weight to each pair of observations
#		"cluster" gives equal weight to each cluster

pvalmax=0

stato=statvc2level(y,cluster,m1,n1,p,weight)
statoriginal=(smat %*% stato[[1]] %*% t(smat))/stato[[2]]

eig=eigen(statoriginal)
maxo=max(eig$values)

if(maxo<=0)
	{
	pvalmax=99
	return(pvalmax)
	}

for(i in 1:npermut)
	{
	clusterperm=relist(sample(1:N),cluster)	
	statp=statvc2level(y,clusterperm,m1,n1,p,weight)
	statperm=(smat %*% statp[[1]] %*% t(smat))/statp[[2]]
	maxperm=max(eigen(statperm)$values)
	pvalmax=pvalmax+(maxperm>maxo)
	}

pvalmax/npermut

}
