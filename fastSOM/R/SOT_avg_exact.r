## Contents: 
#		sot_avg_exact, sot_avg_exact_single, sot_avg_exact_list		  

################################################################################
# calculate the spillover table average exactly, automatic detection whether 'list' (rolling windows) or 'single' use
sot_avg_exact <- function(Sigma,A,ncores=1)
# Sigma: either a covariance matrix or a list thereof
# A: either an array consisting of MA coefficient matrices or a list thereof
# ...: one might specify
#	avg_only: in case of 'single', should the value be only the average or a list of average, minimum, and maximum? 
#	helpers: a list delivered by helpers_avg_exact which takes some time to be calculated and depends only on the number of variables (and therefore might/should be calculated only once)  
# ncores: number of cores, only used in 'list' case
#			missing or 1: no parallelization
#			0: automatic detection of number of cores
#			other integer: use this number 
{
	if (is.list(Sigma)) restructure_list(sot_avg_exact_list(Sigma,A,dim(Sigma[[1]])[1],dim(A[[1]])[3],ncores=ncores))
	else sot_avg_exact_single(Sigma,A,dim(Sigma)[1],dim(A)[3])
}

################################################################################
# calculates the spillover table average exactly, single version  
sot_avg_exact_single <- function(Sigma,A,N=dim(Sigma)[1],H=dim(A)[3],helpers=helpers_avg_exact(N))
{
	# incomplete: handle the case when N is small
	# or do that within the C function (quite unnecessary)
	#if(N<5) return(VDT_DYKW(Sigma=Sigma,A=A,perm=perm)) else require(gtools)
	
	scaling_factor <- 1/sqrt(.Call("fev",Sigma,A,N,H,PACKAGE="fastSOM"))
	res <- .Call("SOT_avg",.Call("scaleSigma",Sigma,scaling_factor,N,PACKAGE="fastSOM"),.Call("scaleA",A,scaling_factor,N,H,PACKAGE="fastSOM"),N,H,helpers$NcK,helpers$cumpos,helpers$gensets-1L,helpers$NminusOne,PACKAGE="fastSOM")
	for (i in 1:3)
	{
		dim(res[[i]]) <- c(N,N)
		dimnames(res[[i]]) <- dimnames(Sigma)
		res[[i]] <- 100*res[[i]]
	}
	res
}

################################################################################
# calculates the spillover table average exactly, list version  
sot_avg_exact_list <- function(Sigma,A,N=dim(Sigma[[1]])[1],H=dim(A[[1]])[3],helpers=helpers_avg_exact(N),ncores=1)
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
			res[[i]] <- sot_avg_exact_single(Sigma[[i]],A[[i]],N,H,helpers)
	}
	else # parallelization
	{
		if (ncores==0)
		{
			ncores <- detectCores() # determine number oc cores
			cat("Number of cores detected:",ncores,"\n")
		}
		
		ncores <- min(ncores,len)
		splitted <- splitIndices(len,ncores) # determine how to distribute the workload 
		cl <- makeCluster(ncores) # create cluster
		clusterEvalQ(cl, library(fastSOM)) # load package fastSOM on every core
		clusterExport(cl,c("Sigma","A","N","H","helpers"),envir=environment()) # send variables to every core
		tmp <- clusterApply(cl,1:ncores,function(ind) sot_avg_exact_list(Sigma[splitted[[ind]]],A[splitted[[ind]]],N,H,helpers,1)) # do parallel jobs
		stopCluster(cl) # close Cluster
		
		for (i in 1:ncores) # putting results together
		{
			res[splitted[[i]]] <- tmp[[i]]
		}
	}
	names(res) <- names(Sigma) 
	res
}
