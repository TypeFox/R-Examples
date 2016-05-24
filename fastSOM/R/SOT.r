## Contents: 
#		sot, sot_single, sot_list

################################################################################
# calculate the spillover table, automatic detection whether 'list' (rolling windows) or 'single' use
sot <- function(Sigma,A,ncores=1,...)
# Sigma: either a covariance matrix or a list thereof
# A: either an array consisting of MA coefficient matrices or a list thereof
# ...: one might specify
#	perm: permutation for which spillover table is to be calculated; if missing, perm will be set to the identity
# ncores: number of cores, only used by 'list' version
#			missing or 1: no parallelization
#			0: automatic detection of number of cores
#			other integer: use this number 
{
	if (is.list(Sigma)) 
		sot_list(Sigma,A,dim(Sigma[[1]])[1],dim(A[[1]])[3],ncores=ncores,...)
	else sot_single(Sigma,A,dim(Sigma)[1],dim(A)[3],...)
}

################################################################################
# calculate the spillover table
sot_single <- function(Sigma,A,N=dim(Sigma)[1],H=dim(A)[3],perm=1:N)
{
	res <- numeric(N*N)
	dim(res) <- c(N,N)
	scaling_factor <- 1/sqrt(.Call("fev",Sigma,A,N,H,PACKAGE="fastSOM"))
	res[] <- 100*.sot_FAST(.Call("scaleSigma",Sigma,scaling_factor,N,PACKAGE="fastSOM"),.Call("scaleA",A,scaling_factor,N,H,PACKAGE="fastSOM"),N,H,perm)
	dimnames(res) <- dimnames(Sigma)
	res
}

################################################################################
# calculates the spillover table, list version  
sot_list <- function(Sigma,A,N=dim(Sigma[[1]])[1],H=dim(A[[1]])[3],perm=1:N,ncores=1)
{
	len <- length(Sigma)
	res <- vector("list",len)
	
	if ( (ncores!=1) && (!require("parallel")) )
	{
		print("Parallelization not possible because package 'parallel' is not installed. Using single core version instead.")
		ncores <- 1
	}
	if (ncores==1) # no parallelization
	{
		for (i in 1:len)
			res[[i]] <- sot_single(Sigma[[i]],A[[i]],N,H,perm)
	}
	else # parallelization
	{
		if (ncores==0)
		{
			ncores <- detectCores() # determine number oc cores
			cat("Number of cores used:",ncores,"\n")
		}
		
		splitted <- splitIndices(len,ncores) # determine how to distribute the workload 
		cl <- makeCluster(ncores) # create cluster
		clusterEvalQ(cl, library(fastSOM)) # load package fastSOM on every core
		clusterExport(cl,c("Sigma","A","N","H","perm"),envir=environment()) # send variables to every core
		tmp <- clusterApply(cl,1:ncores,function(ind) sot_list(Sigma[splitted[[ind]]],A[splitted[[ind]]],N,H,perm,1)) # do parallel jobs
		stopCluster(cl) # close Cluster
		
		for (i in 1:ncores) # putting results together
		{
			res[splitted[[i]]] <- tmp[[i]]
		}
	}
	names(res) <- names(Sigma) 
	res
}

