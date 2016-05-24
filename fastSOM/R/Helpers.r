## Contents: 
# 		.soi_FAST, .soi_FAST_perms, .sot_FAST, .sot_FAST_perms
#		normalize_fev, standardNewN, all_N, helpers_avg_exact, combinations, permutations, 
#		reorganize_list, diagFAST

################################################################################
# handles the use of argument perms in soi_avg_est_single, soi_avg_est_list, sot_avg_est_single, sot_avg_est_list
handle_perms <- function(perms,N)
{
	if (!is.matrix(perms))
	{
		perms <- gen_perms(N,perms)
	}
	else
	if (nrow(perms)!=N)
	{
		if (ncol(perms)==N)
		{
			warning("Dimensions of 'perms' not correct. Will use transposition of perms.")
			perms <- t(perms)
		}
		else 
			stop("perms not suitable!")
	}
	perms
}

################################################################################
# generates nperms random permutations of 1:N
gen_perms <- function(N=19,nperms=1000)
{
	res <-	integer(nperms*N)
	dim(res) <- c(N,nperms)
	for (i in 1:nperms) res[,i] <- sample(N, N, FALSE, NULL)
	res
}

################################################################################
# calculates the spillover index without normalizing, i.e. assuming everything is already normalized
.soi_FAST <- function (Sigma, A, N=dim(Sigma)[2], H=dim(A)[3], perm=1:N) 
{
	.Call("trALsquared",Sigma,A,N,H,perm-1L,PACKAGE="fastSOM")
}

################################################################################
# calculates the spillover index's average over many permutations without normalizing, i.e. assuming everything is already normalized
.soi_FAST_perms <- function (Sigma, A, N=dim(Sigma)[2], H=dim(A)[3], perms, nperms, firstperm=1L) 
{
	tmp <- .Call("trALsquared_perms",Sigma,A,N,H,perms-1L,nperms,firstperm-1L,PACKAGE="fastSOM")
	list(Average=100*(1-tmp[1]/N),
			Min=100*(1-tmp[2]/N),
			Max=100*(1-tmp[3]/N),
			permMin=as.integer(tmp[3+(1:N)]),
			permMax=as.integer(tmp[3+N+(1:N)])
	)
	# using only the average, ignoring max, min, permMax, permMin
}

################################################################################
# calculate the spillover table without normalizing, i.e. assuming everything is already normalized
.sot_FAST <- function(Sigma,A,N=dim(Sigma)[1],H=dim(A)[3],perm=1:N)
{
	.Call("ALsquared",Sigma,A,N,H,perm-1L,PACKAGE="fastSOM")
}

################################################################################
# calculate the spillover table's average over many permutations, without normalizing, i.e. assuming everything is already normalized
.sot_FAST_perms <- function(Sigma,A,N=dim(Sigma)[1],H=dim(A)[3],perms,nperms)
{
	res <- .Call("ALsquared_perms",Sigma,A,N,H,perms-1L,nperms,PACKAGE="fastSOM")
	for (i in 1:3) 
	{
		dim(res[[i]]) <- c(N,N)
		dimnames(res[[i]]) <- dimnames(Sigma)
		res[[i]] <- 100*res[[i]]
	}	
	res
}


################################################################################
# normalize Sigma, A such that forecast error variances are unity
normalize_fev <- function(Sigma,A,N=dim(Sigma)[1],H=dim(A)[3])
{
# rescales Sigma, A to unit forecast error variances 
# Sigma is a covariance matrix, A is an array of MA coefficients
	
	scaling_factor <- 1/sqrt(forecast_error_variances(Sigma,A)) # in order to produce unit forecast error variances 
	
	Sigma[] <- .Call("scaleSigma",Sigma,scaling_factor,N,PACKAGE="fastSOM")
	
	A[] <- .Call("scaleA",A,scaling_factor,N,H,PACKAGE="fastSOM")
	
	list(Sigma=Sigma,A=A)
}

################################################################################
# rule used for 'divide et impera' (calculates the dimension of the 'upper problem')   
standardNewN <- function(n) n %/% 2

################################################################################
# calculates the dimensions that will occur during the 'divide et impera'
all_N <- function(N,Nmin=6,newN=standardNewN)
{
	if (!(N<Nmin))
	{
		N1 <- newN(N)
		N2 <- N-N1
		c(sort(unique(union(Recall(N1,Nmin=Nmin,newN=newN),Recall(N2,Nmin=Nmin,newN=newN)))),N)
	}
	else N
}

################################################################################
# calculate helpers for sot_avg_exact
helpers_avg_exact <- function(N)
{
	NminusOne <- as.integer(.Call("Nminus1",N, PACKAGE="fastSOM"))
	NcK <- as.integer(choose(N-1,0:(N-1)))
	cumpos <- cumsum(c(0L,NcK))
	len <- cumpos[N+1]
	res <- matrix(NA,nrow=len,ncol=N-1)
	res[1,] <- 1:(N-1)
	for (i in 1:(N-1))
	{
		res[(1+cumpos[i+1]):cumpos[i+2],seq_len(i)] <- combinations(N-1,i)
		for (j in 1:NcK[i+1]) res[j+cumpos[i+1],-seq_len(i)] <- (1:(N-1))[-res[j+cumpos[i+1],seq_len(i)]]
	}
	list(NcK=NcK,cumpos=cumpos,gensets=res,NminusOne=NminusOne)
}

################################################################################
# Faster version of combinations, originally defined in gtools, used in helpers_avg_exact
combinations <- function (n, r) 
{
	
	{
		sub <- function(n, r, v) {
			if (r == 1) 
				#.Internal(matrix(v, n, 1, FALSE, NULL, FALSE, FALSE))
				matrix(v,nrow=n,ncol=1)
			else if (r == n) 
				#.Internal(matrix(v, 1, n, FALSE, NULL, FALSE, FALSE))
				matrix(v,nrow=1,ncol=n)
			else rbind(cbind(v[1], Recall(n - 1, r - 1, v[-1])), 
						Recall(n - 1, r, v[-1]))
		}
		sub(n, r, 1:n)
	}
	sub(n, r, 1:n)
}

# Faster version of permutations, originally defined in gtools, used in soi_avg_exact_BruteForce and soi_avg_exact_single
permutations <- function (n) 
{
	sub <- function(n, v) 
	{
		if (n == 1)	
			#.Internal(matrix(v, 1, 1, FALSE, NULL, FALSE, FALSE))
			matrix(v,nrow=1,ncol=1)
		else 
		{
			X <- NULL
			for (i in 1:n) X <- rbind(X, cbind(v[i], Recall(n - 1, v[-i])))
			X
		}
	}
	sub(n, 1:n)
}

################################################################################
# reorganize list of results, used to reorganize results of soi_avg_est_list and soi_avg_exact_list
reorganize_list <- function(data) 
{
	len1 <- length(data)
	len2 <- length(data[[1]])
	res <- vector("list",len2)
	names(res) <- names(data[[1]])
	for ( i in 1:len2)
	{
		#res[[i]] <- vector("list",len1)
		len3 <- length(data[[1]][[i]])
		if (len3>1)
		{
			res[[i]] <- matrix(nrow=len3,ncol=len1)
			colnames(res[[i]]) <- names(data)
			rownames(res[[i]]) <- names(data[[1]][[i]])
			for (j in 1:len1) res[[i]][,j] <- data[[j]][[i]]
		}
		else
		{
			res[[i]] <- vector(typeof(data[[1]][[i]]),len1)
			names(res[[i]]) <- names(data)
			for (j in 1:len1) res[[i]][j] <- data[[j]][[i]]
		}	
	}
	res
}

################################################################################
# reorganize list of results, used to reorganize results of sot_avg_est_list and sot_avg_exact_list
restructure_list <- function(data) 
{
	len1 <- length(data)
	len2 <- length(data[[1]])
	res <- vector("list",len2)
	names(res) <- names(data[[1]])
	for ( i in 1:len2)
	{
		res[[i]] <- vector("list",len1)
		names(res[[i]]) <- names(data)
		for (j in 1:len1) 
		{
			res[[i]][[j]] <- data[[j]][[i]]
		}	
	}	
	res
}

################################################################################
# calculate the diagonal of a matrix, fast version
diagFAST <- function (x = 1) 
{
	m <- dim(x)[1]
	c(x)[1L + 0L:(m - 1L) * (m + 1L)]
}

################################################################################
# calculates the spillover index without normalizing, i.e. assuming everything is already normalized
soi_avg_exact_BruteForce <- function(Sigma,A,B=NULL,N=dim(Sigma)[1],H=dim(A)[3],useB=FALSE,perms=permutations(N))
{
	if (useB)
		.Call("trALplusBLinv_squared_perms",Sigma,A,B,N,H,perms-1L,dim(perms)[2],PACKAGE="fastSOM")
	else	
		.Call("trALsquared_perms",Sigma,A,N,H,perms-1L,dim(perms)[2],0L,PACKAGE="fastSOM")
}


################################################################################
# calculates the spillover index without normalizing, i.e. assuming everything is already normalized
solve_generalized_problem <- function(Sigma,A,N=dim(Sigma)[2],H=dim(A)[3],B=0*A,useB=FALSE,perms,combs,
		Nmin=6,newN=standardNewN,firstlevel=TRUE,ncores=1)
{
# calculates average, variance, minimum und maximum of (the numerator) of the spillover index (without normalizing, i.e. assuming everything is already normalized)  	
	if (N<Nmin)
		return(soi_avg_exact_BruteForce(Sigma,A,B,N,H,useB,perms[[N]]))
	N1 <- newN(N)
	N2 <- N-N1
	tmp <- divide_et_impera(Sigma,A,N,N1,N2,H,B,useB,perms,combs,Nmin,newN,firstlevel,ncores)
	.Call("paste_together",tmp$res1,tmp$res2,N,N1,N2,combs[[N]],dim(combs[[N]])[2],PACKAGE="fastSOM")
}

################################################################################
# calculates the spillover index without normalizing, i.e. assuming everything is already normalized
divide_et_impera <- function(Sigma,A,N,N1,N2,H,B,useB=FALSE,perms,combs,Nmin=6,newN=standardNewN,firstlevel=FALSE,ncores=1)
{
	cur_combs <- combs[[N]]
	ncombs <- dim(cur_combs)[2]
	res1 <- numeric(ncombs*(3+2*N1))
	dim(res1) <- c(3+2*N1,ncombs)
	res2 <- numeric(ncombs*(3+2*N2))
	dim(res2) <- c(3+2*N2,ncombs)
	FirstN1 <- 1:N1
	SecondN2 <- (N1+1):N
	
	parallel <- (firstlevel) && (ncores!=1)
	if ( (parallel) && (!require("parallel")) )
	{
		print("Parallelization not possible because package 'parallel' is not installed. Using single core version instead.")
		ncores <- 1
		parallel <- FALSE
	}
	
	if (!parallel)
	{
		for (i in 1:ncombs)
		{
			M1 <- cur_combs[FirstN1,i]
			M2 <- cur_combs[SecondN2,i]
			Sigma11 <- Sigma[M1,M1]
			Sigma21 <- Sigma[M2,M1]
			Sigma12 <- Sigma[M1,M2]
			
			#tmpM <- solve(Sigma11,Sigma12)#.Call("La_dgesv", Sigma11, Sigma12, 1e-7, PACKAGE = "base") # Sigma11^{-1} %*% Sigma12
			#tmpL <- chol(Sigma11)
			#tmpM <- backsolve(tmpL,forwardsolve(t(tmpL),Sigma12))
			
			.Call("solve_sym",Sigma11,Sigma12,N1,N2,PACKAGE="fastSOM")
			tmpM <- Sigma12
		
			B1 <- B[M1,M1,]
			B2 <- B[M2,M2,]
			if (useB) 
				.Call("array_stuff",B2,-B[M2,M1,],tmpM,N2,N1,H, PACKAGE="fastSOM") # this changes B2[,,h] into B2[,,h] - B[M1,M2,h] %*% tmpM
			.Call("array_stuff",B1,A[M1,M2,],Sigma21,N1,N2,H, PACKAGE="fastSOM") # this changes B1[,,h] into B1[,,h] + A[M1,M2,h] %*% Sigma21
			Sigma22 <- Sigma[M2,M2] 
			.Call("matrix_stuff",Sigma22,-Sigma21,tmpM,N2,N1,N2, PACKAGE="fastSOM") # this changes Sigma22 into Sigma22 - Sigma21 %*% tmpM
			res1[,i] <- solve_generalized_problem(Sigma11,A[M1,M1,],N1,H,B1,TRUE,perms,combs,Nmin,newN,FALSE)
			res2[,i] <- solve_generalized_problem(Sigma22,A[M2,M2,],N2,H,B2,useB,perms,combs,Nmin,newN,FALSE)
		}
	}
	else
	{ # parallel version
		par_problem <- function(cur_combs,Sigma,A,N,N1,N2,H,B,useB,perms,combs,Nmin,newN)
		{
			ncombs <- dim(cur_combs)[2]
			res1 <- numeric(ncombs*(3+2*N1))
			dim(res1) <- c(3+2*N1,ncombs)
			res2 <- numeric(ncombs*(3+2*N2))
			dim(res2) <- c(3+2*N2,ncombs)
			FirstN1 <- 1:N1
			SecondN2 <- (N1+1):N
			for (i in 1:ncombs)
			{
				M1 <- cur_combs[FirstN1,i]
				M2 <- cur_combs[SecondN2,i]
				Sigma11 <- Sigma[M1,M1]
				Sigma21 <- Sigma[M2,M1]
				Sigma12 <- Sigma[M1,M2]
				#tmpM <- solve(Sigma11,Sigma12) 
				#.Call("La_dgesv", Sigma11, Sigma12, 1e-7, PACKAGE = "base") # Sigma11^{-1} %*% Sigma12
				.Call("solve_sym",Sigma11,Sigma12,N1,N2,PACKAGE="fastSOM")
				tmpM <- Sigma12
			
				B1 <- B[M1,M1,]
				B2 <- B[M2,M2,]
				if (useB) 
					.Call("array_stuff",B2,-B[M2,M1,],tmpM,N2,N1,H, PACKAGE="fastSOM") # this changes B2[,,h] into B2[,,h] - B[M1,M2,h] %*% tmpM
				.Call("array_stuff",B1,A[M1,M2,],Sigma21,N1,N2,H, PACKAGE="fastSOM") # this changes B1[,,h] into B1[,,h] + A[M1,M2,h] %*% Sigma21
				Sigma22 <- Sigma[M2,M2] 
				.Call("matrix_stuff",Sigma22,-Sigma21,tmpM,N2,N1,N2, PACKAGE="fastSOM") # this changes Sigma22 into Sigma22 - Sigma21 %*% tmpM
				res1[,i] <- solve_generalized_problem(Sigma11,A[M1,M1,],N1,H,B1,TRUE,perms,combs,Nmin,newN,FALSE)
				res2[,i] <- solve_generalized_problem(Sigma22,A[M2,M2,],N2,H,B2,useB,perms,combs,Nmin,newN,FALSE)
			}
			list(res1=res1,res2=res2)
		}
		if (ncores==0)
		{
			ncores <- detectCores() # determine number of cores
			cat("Number of cores detected:",ncores,"\n")
		}
		
		ncores <- min(ncombs,ncores)
		splitted <- splitIndices(ncombs,ncores) # determine how to distribute the workload 
		cl <- makeCluster(ncores) # create cluster
		clusterEvalQ(cl, library(fastSOM)) # load package fastSOM on every core
		clusterExport(cl,c("cur_combs","Sigma","A","N","N1","N2","H","B","useB","perms","combs","Nmin","newN","par_problem"),envir=environment()) # send variables to every core
		tmp <- clusterApply(cl,1:ncores,function(ind) par_problem(cur_combs[,splitted[[ind]]],Sigma,A,N,N1,N2,H,B,useB,perms,combs,Nmin,newN)) # do parallel jobs
		stopCluster(cl) # close Cluster
		
		for (i in 1:ncores) # putting results together
		{
			res1[,splitted[[i]]] <- tmp[[i]]$res1
			res2[,splitted[[i]]] <- tmp[[i]]$res2
		}
	}
	list(res1=res1,res2=res2)
}


