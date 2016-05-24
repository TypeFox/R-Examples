## Contents: 
#		soi_avg_est, soi_avg_est_single, soi_avg_est_list

################################################################################
# calculate the spillover index average by estimation, automatic detection whether 'list' (rolling windows) or 'single' use
soi_avg_est <- function(Sigma,A,ncores=1,...)
# Sigma: either a covariance matrix or a list thereof
# A: either an array consisting of MA coefficient matrices or a list thereof
# ...: one might specify 
# 	perms: either
#			missing: then 10.000 permutations will be created randomly
#			a matrix: rows must contain permutations of 1,...,N
#			an integer: number of permutations to be created randomly
# ncores: number of cores, only used in 'list' case
#			missing or 1: no parallelization
#			0: automatic detection of number of cores
#			other integer: use this number 
{
	if (is.list(Sigma)) reorganize_list(soi_avg_est_list(Sigma,A,dim(Sigma[[1]])[1],dim(A[[1]])[3],ncores=ncores,...))
	else soi_avg_est_single(Sigma,A,dim(Sigma)[1],dim(A)[3],ncores=ncores,...)
}

################################################################################
# approximates the spillover index average by using 'many' permutations  
soi_avg_est_single <- function(Sigma,A,N=dim(Sigma)[1],H=dim(A)[3],perms,ncores=1)
{
# Sigma is a covariance matrix, A is an array of MA coefficients
# perms, if present, is either
# 	a matrix which contains permutations in its columns; 
#	a number determining the number of permutation to be selected randomly
#	10000, if perms is missing
	if (missing(perms)) 
	{
		perms <- gen_perms(N,10000) 
	} 
	else
	{
		perms <- handle_perms(perms,N)
	}
	nperms <- dim(perms)[2]
	
	tmp <- normalize_fev(Sigma,A,N,H)
	# from here on: parallelization?!
	parallel <- (ncores!=1)
	if ( (parallel) && (!require("parallel")) )
	{
		print("Parallelization not possible because package 'parallel' is not installed. Using single core version instead.")
		ncores <- 1
		parallel <- FALSE
	}
	
	if (parallel)
	{
		if (ncores==0)
		{
			ncores <- detectCores() # determine number of cores
			cat("Number of cores detected:",ncores,"\n")
		}
		
		ncores <- min(ncores,nperms)
		splitted <- splitIndices(nperms,ncores) # determine how to distribute the workload 
		cl <- makeCluster(ncores) # create cluster
		clusterEvalQ(cl, library(fastSOM)) # load package fastSOM on every core
		clusterExport(cl,c("tmp","N","H"),envir=environment()) # send variables to every core ,"perms"
		res <- clusterApply(cl,1:ncores,function(ind) .soi_FAST_perms(tmp$Sigma,tmp$A,N,H,perms=perms,nperms=length(splitted[[ind]]),firstperm=splitted[[ind]][1]) ) # do parallel jobs
		stopCluster(cl) # close Cluster
		
		tmp <- res
		res <- list(Average=0,Min=100,Max=0,permMin=1:N,permMax=1:N)
		for (i in 1:ncores)
		{
			res$Average <- res$Average + tmp[[i]]$Average*length(splitted[[i]])/nperms
			res$Min <- min(res$Min,tmp[[i]]$Min)
			res$Max <- max(res$Max,tmp[[i]]$Max)
			if (res$Min == tmp[[i]]$Min) res$permMin <- tmp[[i]]$permMin 
			if (res$Max == tmp[[i]]$Max) res$permMax <- tmp[[i]]$permMax
		}
		res
	}
	else
	{
		.soi_FAST_perms(tmp$Sigma,tmp$A,N,H,perms,nperms)
	}
}

################################################################################
# approximates the spillover index average by using 'many' permutations  
soi_avg_est_list <- function(Sigma,A,N=dim(Sigma[[1]])[1],H=dim(A[[1]])[3],perms,ncores=1)
{
# Sigma, A are assumed to be lists of covariance matrices and MA coefficients arrays, respectively
# The dimensions of Sigma, A are assumed to be the same for all elements of the list!! 
# perms, if present, is either
# 	a matrix which contains permutations in its rows; 
#	a number determining the number of permutation to be selected randomly
#	10000, if perms is missing
# ncores is the number of cores to be used: 1 (the standard) means no parallelization; if it is set to 0, the number of cores will be detected automatically 
	if (missing(perms)) 
	{
		perms <- gen_perms(N,10000) 
	} 
	else
	{
		perms <- handle_perms(perms,N)
	}
	nperms <- dim(perms)[2]
	
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
			res[[i]] <- soi_avg_est_single(Sigma[[i]],A[[i]],N,H,perms,ncores=1)
	}
	else # parallelization
	{
		if (ncores==0)
		{
			ncores <- detectCores() # determine number oc cores
			cat("Number of cores detected:",ncores,"\n")
		}
		
		ncores <- min(len,ncores)
		splitted <- splitIndices(len,ncores) # determine how to distribute the workload 
		cl <- makeCluster(ncores) # create cluster
		clusterEvalQ(cl, library(fastSOM)) # load package fastSOM on every core
		clusterExport(cl,c("Sigma","A","N","H","perms"),envir=environment()) # send variables to every core
		tmp <- clusterApply(cl,1:ncores,function(ind) soi_avg_est_list(Sigma[splitted[[ind]]],A[splitted[[ind]]],N,H,perms,ncores=1)) # do parallel jobs
		stopCluster(cl) # close Cluster
		
		for (i in 1:ncores) # putting results together
		{
			res[splitted[[i]]] <- tmp[[i]]
		}
	}
	names(res) <- names(Sigma) 
	res
}

