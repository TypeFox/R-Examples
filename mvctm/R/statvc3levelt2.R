statvc3levelt2 <-
function(y,cluster,m1,n1,m2,n2,p,weight){

# y = Matrix Nxp (N obs, p=dimension) of (maybe standardized; signs or ranks) data 
# cluster = List of indices (each element gives the indices of the obs in a level 2 cluster nested within a level 1 cluster). From the "preparecluster" function.
# N = Total number of observations
# p = Dimension
# m1 = Number of level 1 clusters
# n1 = Number of observations in each level 1 cluster (vector m1X1)
# m2 = Number of level 2 clusters within each level 1 cluster (vector m1X1)
# n2 = Number of observations in each level 2 cluster within each level 1 cluster (list of size m1)

mat=matrix(0,p,p)
w=0

for(i in 1:m1)
{
clusteri=cluster[[i]]
m1i=m2[i]
n1i=n2[[i]]
mati=statvc2level(y,clusteri,m1i,n1i,p,weight)
mat=mat+mati[[1]]
w=w+mati[[2]]
}

list(mat,w)
}
