statvc4levelt1 <-
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

if(weight=="observation")
{
for(i in 1:m1)
	{
	clusteri=cluster[[i]]
	m1i=m2[i]
	if(m1i>1)
		{		
		for(j in 1:(m1i-1))
			{
		for(k in (j+1):m1i)
			{
			clusterij=unlist(clusteri[[j]])
			clusterik=unlist(clusteri[[k]])
			nij=length(clusterij)
			nik=length(clusterik)
			wcur=(nij+nik)/2
			w=w+wcur*nij*nik
			yij=y[clusterij,,drop=F]
			yik=y[clusterik,,drop=F]
			sumjk=as.matrix(apply(yij,2,sum)) %*% t(as.matrix(apply(yik,2,sum)))
			mat=mat+wcur*(sumjk+t(sumjk))/2			
			}		
			}
		}
	}
}

else if(weight=="pair")
{
for(i in 1:m1)
	{
	clusteri=cluster[[i]]
	m1i=m2[i]
	if(m1i>1)
		{		
		for(j in 1:(m1i-1))
			{
		for(k in (j+1):m1i)
			{
			clusterij=unlist(clusteri[[j]])
			clusterik=unlist(clusteri[[k]])
			nij=length(clusterij)
			nik=length(clusterik)
			wcur=1
			w=w+wcur*nij*nik
			yij=y[clusterij,,drop=F]
			yik=y[clusterik,,drop=F]
			sumjk=as.matrix(apply(yij,2,sum)) %*% t(as.matrix(apply(yik,2,sum)))
			mat=mat+wcur*(sumjk+t(sumjk))/2				
			}		
			}
		}
	}
}


else if(weight=="cluster")
{
for(i in 1:m1)
	{
	clusteri=cluster[[i]]
	m1i=m2[i]
	if(m1i>1)
		{		
		for(j in 1:(m1i-1))
			{
		for(k in (j+1):m1i)
			{
			clusterij=unlist(clusteri[[j]])
			clusterik=unlist(clusteri[[k]])
			nij=length(clusterij)
			nik=length(clusterik)
			wcur=(nij*nik)
			w=w+wcur*nij*nik
			yij=y[clusterij,,drop=F]
			yik=y[clusterik,,drop=F]
			sumjk=as.matrix(apply(yij,2,sum)) %*% t(as.matrix(apply(yik,2,sum)))
			mat=mat+wcur*(sumjk+t(sumjk))/2				
			}		
			}
		}
	}
}


list(mat,w)
}
