statvc2level <-
function(y,cluster,m1,n1,p,weight){

# cluster = list of indices (each element gives the indices of the obs in a level 1 cluster)
# m1 = number of level 1 cluster (=number of elements in cluster)
# n1 = number of observations in each level 1 cluster (vector m1X1)

mat=matrix(0,p,p)
w=0

if(weight=="observation")
	{
		for(i in 1:m1)
		{
		if(n1[i]>1)
			{
			wcur=1/(n1[i]-1)
			ycl=y[cluster[[i]],,drop=F]
			w=w+wcur*(n1[i]-1)*(n1[i])/2
			matbar=t(as.matrix(apply(ycl,2,mean)))
			mat2=t(ycl) %*% ycl
			mat=mat+wcur*(n1[i]^2*(t(matbar) %*% matbar) - mat2)/2
	}   
		}
			}


else if(weight=="pair")
	{
		for(i in 1:m1)
		{
		if(n1[i]>1)
			{
			wcur=1
			ycl=y[cluster[[i]],,drop=F]
			w=w+wcur*(n1[i]-1)*(n1[i])/2
			matbar=t(as.matrix(apply(ycl,2,mean)))
			mat2=t(ycl) %*% ycl
			mat=mat+wcur*(n1[i]^2*(t(matbar) %*% matbar) - mat2)/2
			}   
		}
	}


else if(weight=="cluster")
	{
		for(i in 1:m1)
		{
		if(n1[i]>1)
			{
			wcur=1/(n1[i]*(n1[i]-1))
			ycl=y[cluster[[i]],,drop=F]
			w=w+wcur*(n1[i]-1)*(n1[i])/2
			matbar=t(as.matrix(apply(ycl,2,mean)))
			mat2=t(ycl) %*% ycl
			mat=mat+wcur*(n1[i]^2*(t(matbar) %*% matbar) - mat2)/2
	}   
		}
			}


list(mat,w)
}
