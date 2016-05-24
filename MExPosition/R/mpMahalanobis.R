## computes Mahalanobis distance
mpMahalanobis <- function(data, row.design)
{	num.groups = dim(row.design)[2]
	num.obs = dim(data)[1]

	Xbar <- matrix(0,dim(row.design)[2],dim(data)[2])	
	for(i in 1:dim(row.design)[2])
    { Xbar[i,] <- colMeans(data[c(which(row.design[,i]==1)),])	
   	}

   	#repeat means per row group	
	Dw <- data - Xbar[rep(1:dim(row.design)[2],c(colSums(row.design))),]

	sigma <- (1/(num.obs - dim(row.design)[2])) * (t(Dw) %*% Dw)
	
	S <- Xbar %*% solve(sigma) %*% t(Xbar)
	s <- diag(S)
	D <- repmat(s,1,num.groups) + repmat(t(s),num.groups,1) - (2*S)
	
	return(D)
}