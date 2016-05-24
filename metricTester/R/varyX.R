#' Calculate alpha or beta metrics across a set of parameters
#'
#' Takes a specified set of parameters, and holds those constant while varying one of
#' the parameters to see variance in community structure metrics.
#'
#' @param alpha Whether to calculate alpha or beta phylogenetic community structure
#' metrics. Default is TRUE. Set to FALSE for beta metrics.
#' @param tree.size Either a single number, or if the aim is to vary tree size while
#' holding other parameters constant, then a list of numbers of species desired in each
#' total tree. If provided as a vector of numbers, will be coerced to a list.
#' @param richness Either the number of species to be placed in each plot, or if the aim
#' is to vary richness while holding other parameters constant, then a list of vectors of
#' number of species to be placed in each plot. See examples.
#' @param delta Either a value for the delta transformation (Pagel 1999), or if the aim is
#' to vary tree shape while holding other parameters constant, then a list of numbers of
#' delta transformations. If provided as a vector of numbers, will be coerced to a list.
#' Values greater than 1 push the branching events towards the root, while values less
#' than 1 push the branching events closer to the tips. See details for particularly low
#' delta values.
#' @param abundances Either a vector of abundances, or if the aim is to vary the abundance
#' distribution function, then a list of vectors of abundances. See examples.
#' @param beta.iterations Because the type of beta-level phylogenetic community structure
#' metrics used here return a single value per community data matrix, it is not possible
#' to look for inter-metric correlations with only a single matrix and tree. To deal with
#' this, the same tree can be used with different community data matrices. This argument
#' specifies the number of matrices to be used per tree. Not needed if alpha=TRUE.
#' @param iterations How many times to simulate the given set of parameters. For instance,
#' with a single tree size, richness.vector, delta, and two sets of abundances, and 10
#' iterations, the result will be a list with 10 elements. Each of those 10 elements will
#' be a list of two elements, each of which will be the calculated metrics for a given
#' set of parameters (one for each abundance vector).
#' @param cores This function can run in parallel. In order to do so, the user must
#' specify the desired number of cores to utilize. The default is "seq", which runs the
#' calculations sequentially. The iteration aspect
#' of the function is parallelized, so for efficiency purposes, it is best to run this
#' over numerous iterations rather than repeating the same parameter numerous times (e.g.,
#' rather than setting deltas to rep(1, 10), set delta to 1 and iterations to 10). 
#' 
#' @details If given a small value, e.g. 0.1, the delta parameter
#' (tree shape) can occasionally result in oddly formatted trees that would cause errors.
#' To deal with this, there is an internal check that will recreate a new tree and
#' re-scale it with the desired delta. This has not been tested at delta < 0.1, and is
#' currently programmed with a while loop. Care should be taken not to get R stuck in an
#' indefinite loop at delta values even lower than 0.1
#'
#' @return A list of lists of data frames.
#'
#' @export
#'
#' @importFrom foreach foreach %dopar% %do%
#' @importFrom doParallel registerDoParallel
#'
#' @references Miller, E. T., D. R. Farine, and C. H. Trisos. 2015. Phylogenetic community
#' structure metrics and null models: a review with new methods and software.
#' bioRxiv 025726.
#'
#' @examples
#' #example of how to vary tree size
#' #not run
#' #system.time(temp <- varyX(alpha=TRUE, tree.size=c(59, 100),
#'	#richness=40:59, delta=1,
#'	#abundances=round(rlnorm(5000, meanlog=2, sdlog=1)) + 1, iterations=2))
#'
#' #example of how to vary richness
#' #not run
#' #system.time(temp <- varyX(alpha=TRUE, tree.size=59,
#'	#richness=list(30:39, 40:49), delta=1,
#'	#abundances=round(rlnorm(5000, meanlog=2, sdlog=1)) + 1, iterations=2))
#'
#' #example of how to vary tree shape
#' #not run
#' #system.time(temp <- varyX(alpha=TRUE, tree.size=59,
#'	#richness=40:59, delta=c(0.1,10),
#'	#abundances=round(rlnorm(5000, meanlog=2, sdlog=1)) + 1, iterations=2))
#'
#' #example of how to vary abundance
#' #not run
#' #inputAbunds <- list(rep(2,5000), (round(rnorm(5000, 50, 15)) + 1),
#'	#(round(rlnorm(5000, meanlog=2, sdlog=1)) + 1))
#'
#' #system.time(temp <- varyX(alpha=TRUE, tree.size=59,
#'	#richness=40:59, delta=1, abundances=inputAbunds, iterations=2))

varyX <-  function(alpha=TRUE, tree.size, richness, delta, abundances,
	beta.iterations, iterations, cores="seq")
{
	#if the inputs match the expectations for varyTreeSize
	if(class(tree.size) == "list" | length(tree.size) > 1
		& class(richness) != "list" & class(delta) != "list" & length(delta) == 1
		& class(abundances) != "list")
	{
		print("Varying tree sizes and holding other parameters constant")
		results <- varyTreeSize(alpha=alpha, tree.sizes=tree.size,
			richness.vector=richness, delta=delta, abundances=abundances,
			beta.iterations=beta.iterations, iterations=iterations, cores=cores)
	}
	#if the inputs match the expectations for varyRichness
	else if(class(tree.size) != "list" & length(tree.size) == 1
		& class(richness) == "list" & class(delta) != "list" & length(delta) == 1
		& class(abundances) != "list")
	{
		print("Varying richnesses and holding other parameters constant")
		results <- varyRichness(alpha=alpha, tree.size=tree.size,
			richness.vectors=richness, delta=delta, abundances=abundances,
			beta.iterations=beta.iterations, iterations=iterations, cores=cores)
	}
	#if the inputs match the expectations for varyTreeShape
	else if(class(delta) == "list" | length(delta) > 1
		& class(richness) != "list" & class(tree.size) != "list" & length(tree.size) == 1
		& class(abundances) != "list")
	{
		print("Varying tree shapes and holding other parameters constant")
		results <- varyTreeShape(alpha=alpha, tree.size=tree.size,
			richness.vector=richness, deltas=delta, abundances=abundances,
			beta.iterations=beta.iterations, iterations=iterations, cores=cores)
	}
	#if the inputs match the expectations for varyAbundance
	else if(class(tree.size) != "list" & length(tree.size) == 1
		& class(richness) != "list" & class(delta) != "list" & length(delta) == 1
		& class(abundances) == "list")
	{
		print("Varying abundance distributions and holding other parameters constant")
		results <- varyAbundance(alpha=alpha, tree.size=tree.size,
			richness.vector=richness, delta=delta, multi.abundances=abundances,
			beta.iterations=beta.iterations, iterations=iterations, cores=cores)

	}
	else
	{
		stop("Carefully check the classes and formatting of your inputs")
	}
	results
}

varyTreeSize <- function(alpha=TRUE, tree.sizes, richness.vector, delta, abundances,
	beta.iterations, iterations, cores)
{
	#this function is not exported. it is called by varyX.
	#if the input is not a list, coerce to list
	if(!is.list(tree.sizes))
	{
		tree.sizes <- as.list(tree.sizes)
	}

	#if alpha is true, run the alpha diversity phylo metric pipeline	
	if(alpha==TRUE)
	{
		#if cores is set to seq, do not register parallel cores
		if(cores == "seq")
		{
			#warn that the analysis is being run sequentially
			warning("Not running analysis in parallel. See 'cores' argument.",
				call.=FALSE)
			
			results <- foreach(i = 1:iterations) %do%
			{
				lapply(seq_along(tree.sizes), function(x)
					alphaMetricSims(tree.size=tree.sizes[[x]], 
					richness.vector=richness.vector, delta=delta, abundances))
			}
		}
		#if cores is not set to seq, spawn the number of cores requested. note that there
		#are currently no checks to ensure the input is a valid number of cores or "seq"
		if(cores != "seq")
		{
			registerDoParallel(cores)
		
			results <- foreach(i = 1:iterations) %dopar%
			{
				lapply(seq_along(tree.sizes), function(x)
					alphaMetricSims(tree.size=tree.sizes[[x]], 
					richness.vector=richness.vector, delta=delta, abundances))
			}	
			registerDoSEQ()
		}
		for(i in 1:length(results))
		{
			names(results[[i]]) <- paste("tree.size", 1:length(tree.sizes), sep="")
		}
	}

	#if alpha is false, run the beta diversity metrics pipeline
	else if(alpha==FALSE)
	{
		if(cores == "seq")
		{
			#warn that the analysis is being run sequentially
			warning("Not running analysis in parallel. See 'cores' argument.",
				call.=FALSE)
			
			results <- foreach(i = 1:iterations) %do%
			{
				lapply(seq_along(tree.sizes), function(x)
					betaMetricSims(tree.size=tree.sizes[[x]], 
					richness.vector=richness.vector, delta=delta, abundances,
					beta.iterations))
			}
		}

		#if cores is not set to seq, spawn the number of cores requested. note that there
		#are currently no checks to ensure the input is a valid number of cores or "seq"
		if(cores != "seq")
		{
			registerDoParallel(cores)	

			results <- foreach(i = 1:iterations) %dopar%
			{
				lapply(seq_along(tree.sizes), function(x)
					betaMetricSims(tree.size=tree.sizes[[x]], 
					richness.vector=richness.vector, delta=delta, abundances,
					beta.iterations))
			}
			registerDoSEQ()
		}

		for(i in 1:length(results))
		{
			names(results[[i]]) <- paste("tree.size", 1:length(tree.sizes), sep="")
		}
	}
	
	names(results) <- paste("iteration", 1:iterations, sep="")
	
	results
}

varyRichness <- function(alpha=TRUE, tree.size, richness.vectors, delta, abundances,
	beta.iterations, iterations, cores)
{
	#this function is not exported. it is called by varyX.
	#if the input is not a list, coerce to list
	if(!is.list(richness.vectors))
	{
		richness.vectors <- as.list(richness.vectors)
	}

	#if alpha is true, run the alpha diversity phylo metric pipeline	
	if(alpha==TRUE)
	{
		#if cores is set to seq, do not register parallel cores
		if(cores == "seq")
		{
			#warn that the analysis is being run sequentially
			warning("Not running analysis in parallel. See 'cores' argument.",
				call.=FALSE)

			results <- foreach(i = 1:iterations) %do%
			{
				lapply(seq_along(richness.vectors), function(x)
					alphaMetricSims(tree.size=tree.size, 
					richness.vector=richness.vectors[[x]], delta=delta, abundances))
			}
		}
		
		if(cores != "seq")
		{
			registerDoParallel(cores)

			results <- foreach(i = 1:iterations) %dopar%
			{
				lapply(seq_along(richness.vectors), function(x)
					alphaMetricSims(tree.size=tree.size, 
					richness.vector=richness.vectors[[x]], delta=delta, abundances))
			}

			registerDoSEQ()
		}

		for(i in 1:length(results))
		{
			names(results[[i]]) <- paste("richness.vector", 1:length(richness.vectors),
				sep="")
		}
	}

	#if alpha is false, run the beta diversity metrics pipeline
	else if(alpha==FALSE)
	{	
		#if cores is set to seq, do not register parallel cores
		if(cores == "seq")
		{
			#warn that the analysis is being run sequentially
			warning("Not running analysis in parallel. See 'cores' argument.",
				call.=FALSE)

			results <- foreach(i = 1:iterations) %do%
			{
				lapply(seq_along(richness.vectors), function(x)
					betaMetricSims(tree.size=tree.size, 
					richness.vector=richness.vectors[[x]], delta=delta, abundances,
					beta.iterations))
			}
		}

		if(cores != "seq")
		{
			registerDoParallel(cores)
			
			results <- foreach(i = 1:iterations) %dopar%
			{
				lapply(seq_along(richness.vectors), function(x)
					betaMetricSims(tree.size=tree.size, 
					richness.vector=richness.vectors[[x]], delta=delta, abundances,
					beta.iterations))
			}
			
			registerDoSEQ()
		}

		for(i in 1:length(results))
		{
			names(results[[i]]) <- paste("richness.vector", 1:length(richness.vectors),
				sep="")
		}
	}
	
	names(results) <- paste("iteration", 1:iterations, sep="")
	
	results
}

varyTreeShape <- function(alpha=TRUE, tree.size, richness.vector, deltas, abundances,
	beta.iterations, iterations, cores)
{
	#this function is not exported. it is called by varyX.
	#if the input is not a list, coerce to list
	if(!is.list(deltas))
	{
		deltas <- as.list(deltas)
	}

	#if alpha is true, run the alpha diversity phylo metric pipeline	
	if(alpha==TRUE)
	{
		#if cores is set to seq, do not register parallel cores
		if(cores == "seq")
		{
			#warn that the analysis is being run sequentially
			warning("Not running analysis in parallel. See 'cores' argument.",
				call.=FALSE)

			results <- foreach(i = 1:iterations) %do%
			{
				lapply(seq_along(deltas), function(x)
					alphaMetricSims(tree.size=tree.size, 
					richness.vector=richness.vector, delta=deltas[[x]], abundances))
			}
		}

		if(cores != "seq")
		{
			registerDoParallel(cores)
			
			results <- foreach(i = 1:iterations) %dopar%
			{
				lapply(seq_along(deltas), function(x)
					alphaMetricSims(tree.size=tree.size, 
					richness.vector=richness.vector, delta=deltas[[x]], abundances))
			}

			registerDoSEQ()
		}

		for(i in 1:length(results))
		{
			names(results[[i]]) <- paste("tree.shape", 1:length(deltas), sep="")
		}
	}

	#if alpha is false, run the beta diversity metrics pipeline
	else if(alpha==FALSE)
	{
		#if cores is set to seq, do not register parallel cores
		if(cores == "seq")
		{
			#warn that the analysis is being run sequentially
			warning("Not running analysis in parallel. See 'cores' argument.",
				call.=FALSE)

			results <- foreach(i = 1:iterations) %do%
			{
				lapply(seq_along(deltas), function(x)
					betaMetricSims(tree.size=tree.size, 
					richness.vector=richness.vector, delta=deltas[[x]], abundances,
					beta.iterations))
			}
		}

		if(cores != "seq")
		{
			registerDoParallel(cores)

			results <- foreach(i = 1:iterations) %dopar%
			{
				lapply(seq_along(deltas), function(x)
					betaMetricSims(tree.size=tree.size, 
					richness.vector=richness.vector, delta=deltas[[x]], abundances,
					beta.iterations))
			}

			registerDoSEQ()
		}

		for(i in 1:length(results))
		{
			names(results[[i]]) <- paste("tree.shape", 1:length(deltas), sep="")
		}
	}
	
	names(results) <- paste("iteration", 1:iterations, sep="")
	
	results
}

varyAbundance <- function(alpha=TRUE, tree.size, richness.vector, delta, multi.abundances,
	beta.iterations, iterations, cores)
{
	#this function is not exported. it is called by varyX.
	#if the input is not a list, coerce to list
	if(!is.list(multi.abundances))
	{
		multi.abundances <- as.list(multi.abundances)
	}

	#if alpha is true, run the alpha diversity phylo metric pipeline	
	if(alpha==TRUE)
	{
		#if cores is set to seq, do not register parallel cores
		if(cores == "seq")
		{
			#warn that the analysis is being run sequentially
			warning("Not running analysis in parallel. See 'cores' argument.",
				call.=FALSE)

			results <- foreach(i = 1:iterations) %do%
			{
				lapply(seq_along(multi.abundances), function(x)
					alphaMetricSims(tree.size=tree.size, 
					richness.vector=richness.vector, delta=delta, multi.abundances[[x]]))
			}
		}
		
		#if cores is not set to seq, spawn the number of cores requested. note that there
		#are currently no checks to ensure the input is a valid number of cores or "seq"
		if(cores != "seq")
		{
			registerDoParallel(cores)	

			results <- foreach(i = 1:iterations) %dopar%
			{
				lapply(seq_along(multi.abundances), function(x)
					alphaMetricSims(tree.size=tree.size, 
					richness.vector=richness.vector, delta=delta, multi.abundances[[x]]))
			}
			
			registerDoSEQ()
		}

		for(i in 1:length(results))
		{
			names(results[[i]]) <- paste("abundance", 1:length(multi.abundances), sep="")
		}
	}

	#if alpha is false, run the beta diversity metrics pipeline
	else if(alpha==FALSE)
	{
		#if cores is set to seq, do not register parallel cores
		if(cores == "seq")
		{
			#warn that the analysis is being run sequentially
			warning("Not running analysis in parallel. See 'cores' argument.",
				call.=FALSE)

			results <- foreach(i = 1:iterations) %do%
			{
				lapply(seq_along(multi.abundances), function(x)
					betaMetricSims(tree.size=tree.size, 
					richness.vector=richness.vector, delta=delta, multi.abundances[[x]],
					beta.iterations))
			}
		}

		#if cores is not set to seq, spawn the number of cores requested. note that there
		#are currently no checks to ensure the input is a valid number of cores or "seq"
		if(cores != "seq")
		{
			registerDoParallel(cores)	

			results <- foreach(i = 1:iterations) %dopar%
			{
				lapply(seq_along(multi.abundances), function(x)
					betaMetricSims(tree.size=tree.size, 
					richness.vector=richness.vector, delta=delta, multi.abundances[[x]],
					beta.iterations))
			}
			
			registerDoSEQ()
		}

		for(i in 1:length(results))
		{
			names(results[[i]]) <- paste("abundance", 1:length(multi.abundances), sep="")
		}
	}
	
	names(results) <- paste("iteration", 1:iterations, sep="")
	
	results
}
