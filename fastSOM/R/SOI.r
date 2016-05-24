## Contents: 
#		soi, soi_single, soi_list, soi_from_sot  

################################################################################
# calculate the spillover index, automatic detection whether 'list' (rolling windows) or 'single' use
soi <- function(Sigma,A,ncores=1,...)
{
	# Sigma: either a covariance matrix or a list thereof
	# A: either an array consisting of MA coefficient matrices or a list thereof
	# ...: one might specify
	#	perm: permutation for which spillover index is to be calculated; if missing, perm will be set to the identity
	# ncores: number of cores, only used in 'list' case
	#			missing or 1: no parallelization
	#			0: automatic detection of number of cores
	#			other integer: use this number 
	if (is.list(Sigma)) soi_list(Sigma,A,dim(Sigma[[1]])[1],dim(A[[1]])[3],ncores=ncores,...) 
	else soi_single(Sigma,A,dim(Sigma)[1],dim(A)[3],...) 
}

################################################################################
# calculates the spillover index, intended for single use
soi_single <- function(Sigma,A,N=dim(Sigma)[1],H=dim(A)[3],perm=1:N)
{
	# first: normalize to unit forecast error variances 
	scaling_factor <- 1/sqrt(forecast_error_variances(Sigma,A))
	100*(1-.soi_FAST(.Call("scaleSigma",Sigma,scaling_factor,N,PACKAGE="fastSOM"),.Call("scaleA",A,scaling_factor,N,H,PACKAGE="fastSOM"),N,H,perm)/N)
}

################################################################################
# calculates the spillover index for lists
soi_list <- function(Sigma,A,N=dim(Sigma[[1]])[1],H=dim(A[[1]])[3],perm=1:N,ncores=1)
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
			res[[i]] <- soi_single(Sigma[[i]],A[[i]],N,H,perm)
	}
	else # parallelization
	{
		if (ncores==0)
		{
			ncores <- detectCores() # determine number of cores
			cat("Number of cores used:",ncores,"\n")
		}
		
		splitted <- splitIndices(len,ncores) # determine how to distribute the workload 
		cl <- makeCluster(ncores) # create cluster
		clusterEvalQ(cl, library(fastSOM)) # load package fastSOM on every core
		clusterExport(cl,c("Sigma","A","N","H","perm"),envir=environment()) # send variables to every core
		tmp <- clusterApply(cl,1:ncores,function(ind) soi_list(Sigma[splitted[[ind]]],A[splitted[[ind]]],N,H,perm,1)) # do parallel jobs
		stopCluster(cl) # close Cluster
		
		for (i in 1:ncores) # putting results together
		{
			res[splitted[[i]]] <- tmp[[i]]
		}
	}
	names(res) <- names(Sigma) 
	res
}

################################################################################
# calculate the spillover index from a spillover table, automatic detection whether 'list' (rolling windows) or 'single' use
soi_from_sot <- function(input_table)
{
	# input_table: either a spillover table or a list thereof
	if (is.list(input_table)) lapply(input_table,function(input) 100-mean(diagFAST(input))) 
	else 100-mean(diagFAST(input_table)) 
}
