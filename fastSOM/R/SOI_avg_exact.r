## Contents: 
#		soi_avg_exact, soi_avg_exact_single, !!soi_avg_exact_list!!  

################################################################################
# calculate the spillover index average by estimation, automatic detection whether 'list' (rolling windows) or 'single' use
soi_avg_exact <- function(Sigma,A,ncores=1)
# Sigma: either a covariance matrix or a list thereof
# A: either an array consisting of MA coefficient matrices or a list thereof
# ...: one might specify 
# ncores: number of cores, only used in 'list' case
#			missing or 1: no parallelization
#			0: automatic detection of number of cores
#			other integer: use this number 
{
	if (is.list(Sigma)) reorganize_list(soi_avg_exact_list(Sigma,A,dim(Sigma[[1]])[1],dim(A[[1]])[3],ncores=ncores))
	else soi_avg_exact_single(Sigma,A,dim(Sigma)[1],dim(A)[3],ncores=ncores)
}

################################################################################
# calculates average, minimum, and maximum of the spillover index, single version, possibly uses parallelization (by supplying ncores!=1)  
soi_avg_exact_single <- function(Sigma,A,N=dim(Sigma)[2],H=dim(A)[3],Nmin=6,NewN=standardNewN,...)
{
	#tmp <- fastSOM:::soi_avg_exact_single(Sigma,A)
	Nmin <- max(Nmin,4)
	vecN <- all_N(N,Nmin,standardNewN) # vector containing the N's occuring during the divide-and-conquer
	perms <- list() 
	combs <- list()
	smallN <- (vecN<Nmin) # indicator for N's that will be solved by brute force
	for (i in vecN[smallN]) # create permutations for N's later on solved by brute force
		perms[[i]] <- t(permutations(i))
	if (any(!smallN)) for (i in vecN[!smallN]) # create combinations for N's later on solved by divide-and-conquer
		{
			combs[[i]] <- t(combinations(i,NewN(i)))
			combs[[i]] <- rbind(combs[[i]],apply(combs[[i]],2,function(x) (1:i)[-x]))
		}	
	
	# scaling
	tmp <- normalize_fev(Sigma,A,N,H)
	# solve generalized problem
	res <- solve_generalized_problem(tmp$Sigma,tmp$A,N,H,B=0*tmp$A,useB=FALSE,perms,combs,Nmin,NewN,...)
	
	# calculate return values
	list(Average=100*(1-res[1]/N),Min=100*(1-res[2]/N),Max=100*(1-res[3]/N),
			permMin=as.integer(res[3+(1:N)]),permMax=as.integer(res[3+N+(1:N)]))
}

################################################################################
# approximates the spillover index average by using 'many' permutations  
soi_avg_exact_list <- function(Sigma,A,N=dim(Sigma[[1]])[1],H=dim(A[[1]])[3],ncores=1,...)
{
# Sigma, A are assumed to be lists of covariance matrices and MA coefficients arrays, respectively
# The dimensions of Sigma, A are assumed to be the same for all elements of the list!! 
# ncores is the number of cores to be used: 1 (the standarad) means no parallelization; if it is set to 0, the number of cores will be detected automatically 
	
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
			res[[i]] <- soi_avg_exact_single(Sigma[[i]],A[[i]],N,H,ncores=1,...)
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
		tmp <- clusterApply(cl,1:ncores,function(ind) soi_avg_exact_list(Sigma[splitted[[ind]]],A[splitted[[ind]]],N,H,ncores=1,...)) # do parallel jobs
		stopCluster(cl) # close Cluster
		
		for (i in 1:ncores) # putting results together
		{
			res[splitted[[i]]] <- tmp[[i]]
		}
	}
	names(res) <- names(Sigma) 
	res
}

