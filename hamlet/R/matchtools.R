#
# Generalized tools for bipartite and non-bipartite matching
#


# A generic framework for the genetic algorithm for performing multigroup non-bipartite matching (local optimum)
match.ga <- function(
	# Distance/dissimilarity matrix d
	d,
	# Subgroup size; g=2 is pairs, g=3 is triplets, etc
	g,
	# If user wants to provide existing population matrix, pops-parameter can be used
	pops,
	# Number of generations
	generations = 100,
	# Population size, i.e. number of solutions simulated per each generation
	popsize = 100,
	# Number of mutations per generation
	nmutate = 100,
	# Number of solution deaths (and thus breeding) per generation
	ndeath = 30,
	# Type, do we desire minimization ("min") or maximization ("max")
	type="min",
	# Mutation function, how point mutations happen for a solution vector 'x'	
	# Should return an eligible solution vector
	mutate = hamlet:::.ga.mutate,
	# How breeding is done when two solution vectors 'x' and 'y' are combined to produce offspring
	# Should return an eligible solution vector	
	breed = hamlet:::.ga.breed
	,
	# How to compute weighting for probabilities to survive / breed / get mutations
	weight = hamlet:::.ga.weight
	,	
	# The way the fitness value of a solution is computed - most likely the value of the function to be optimized
	# Input 'x' is a solution vector of submatches per each individual
	# Should return a single value that describes fitness of a solution (higher -> better)
	fitness = hamlet:::.ga.fitness
	,
	# Step function that produces the next generation of solutions
	# pops: Matrix of previous generation solutions (rows = individuals, cols = different matching solution vectors)
	# fitnesses: Vector of fitness values for the solutions (length = number of columns in 'pops')
	step = hamlet:::.ga.step
	,
	# Generate random starting conditions; a matrix with nrow(d), ncol(d) amount of rows, and 'popsize' amount of columns
	# Thus, columns depict the solution vectors living in the genetic algorithm
	# Rows are the actual individuals in the study
	initialize = hamlet:::.ga.init
	,
	# If plots should be given during algorithm
	progplot = T,
	# If a generation-plot should be done to illustrate advancement of the algorithm
	plot = T,
	# Level of verbosity; <0 means completely silent, >=0 means normal verbosity, >=1 means also additional info
	verb = 0,
	# After how many generations a new progress plot is done
	progress = 500,
	...
){
	# Make sure d is a matrix
	d <- as.matrix(d)
	# Make sure g is an integer
	g <- round(g, 0)
	# Handling invalid input 
	if(missing(g)) stop("Parameter 'g' missing, indicate subgroup size with an integer")
	if(!nrow(d) %% g == 0 | !ncol(d) %% g == 0) stop("Number of rows and columns in 'd' should be dividable by 'g'")
	if(!nrow(d) == ncol(d)) stop("Distance matrix 'd' should be a square matrix")

	sols <- list()
	best <- ifelse(type=="min", Inf, -Inf) # Cost of best found solution so far
	bestsol <- NA # Best found solution so far
	quants <- matrix(nrow=5, ncol=generations) # Matrix for storing and plotting solution quantiles per time
	rownames(quants) <- c("0%", "25", "50%", "75%", "100%")
	# Initial values (initial population matrix can also be supplied by the user with param 'pops')	
	if(missing(pops)) pops <- initialize(popsize = popsize, d = d, g = g)
	# Compute fitnesses and weights
	fitnesses <- apply(pops, MARGIN=2, FUN=function(z) fitness(z, d=d))
	weights <- weight(fitnesses)
	# Start iterating through the genetic algorithm
	for(i in 1:generations){
		sols <- list(pops = pops, fitnesses = fitnesses, weights = weights)
		new <- step(pops = pops, fitnesses = fitnesses, weights = weights, nmutate = nmutate, ndeath = ndeath, mutate = mutate, breed = breed, fitness = fitness, weight = weight, d = d)
		pops <- new$pops
		fitnesses <- new$fitnesses
		weights <- new$weights
		quants[,i] <- quantile(fitnesses)
		# New best found solution
		if(ifelse(type=="min", quants[1,i]<best, quants[1,i]>best)){
			best <- ifelse(type=="min", quants[1,i], quants[5,i])
			bestsol <- pops[,which(fitnesses==best)[1]]
		}
		# Print progress information
		if(i %% progress == 0){ 
			print(paste("Generation", i, "of", generations)) 			
			print("Current quantiles:") 
			print(quants[,i])
			print("Best found solution vector:")
			print(bestsol)
			print("Best found solution cost:")
			print(best)
			# Plot progress quantiles
			if(progplot){
				plot(quants[3,1:i], type="l", ylim=extendrange(quants[,1:i]), xlab="Generation", ylab="Fitness")
				points(quants[1,1:i], type="l", col="red")
				points(quants[5,1:i], type="l", col="red")
				points(quants[2,1:i], type="l", col="orange", lty="dashed")
				points(quants[4,1:i], type="l", col="orange", lty="dashed")
				legend("topleft", horiz=T, bty="n", col=c("red", "orange", "black"), lty=c("solid", "dashed", "solid"), legend=c("Min/Max", "2/4 Quant.", "Median"), cex=0.7)
			}
		}
	}
	# Last iteration collecting from the results in the for-loop
	sols <- list(pops = pops, fitnesses = fitnesses, weights = weights)
	
	# Plot the final progress plot
	if(plot){
		plot(quants[3,], type="l", ylim=extendrange(quants[,1:i]), xlab="Generation", ylab="Fitness")
		points(quants[1,], type="l", col="red")
		points(quants[5,], type="l", col="red")
		points(quants[2,], type="l", col="orange", lty="dashed")
		points(quants[4,], type="l", col="orange", lty="dashed")
		legend("topleft", horiz=T, bty="n", col=c("red", "orange", "black"), lty=c("solid", "dashed", "solid"), legend=c("Min/Max", "2/4 Quant.", "Median"), cex=0.7)
	}
	if(verb>=0){
		print("Best found solution vector:")
		print(bestsol)
		print("Best found solution cost:")
		print(best)
	}

	list(sols, bestsol, best)
}

# Branch & Bound algorithm for multigroup non-bipartite optimal matching (global optimum)
match.bb <- function(
	# Symmetric distance/dissimilarity matrix d
	d,
	# Group size, how many individuals belong to each sub-match
	g = 2,
	# Method for presorting the distance matrix; values should be lower near diagonal
	# Default is complete hierarchic clustering, possible values: "complete", "ward", "single", "average", "mcquitty", "median", "centroid"
	presort = "complete",
	# After how many branching operations we will inform user about the current progress
	progress = 100000,
	# Starting solution; a value of a solution that is known to exist, the lower the better (more bounding)
	# if none are known, this should be set to Inf, however even a naive preliminary solution may be nice
	bestknown = Inf,
	# How many branching operations are allowed before stopping to possibly local optimum
	maxbranches = Inf,
	# Level of verbosity, real number; lower means less, higher means more
	verb = 0
){
	# Test feasibility of 'd' for the matching task
	if(nrow(d)<2 | ncol(d)<2)
		stop("Distance/dissimilarity matrix 'd' should be at least 2x2")
	if(!dim(d)[1]%%g ==0 | !dim(d)[2]%%g==0)
		stop("Dimensions of the distance/dissimilarity matrix 'd' should be dividable by 'g', consider adding dummy individuals")
	if(!dim(d)[1]==dim(d)[2])
		stop("Distance/dissimilarity should be a square matrix")

	# Variables etc
	# Setting diagonal to NA
	diag(d) <- NA
	nr = nrow(d)
	nc = ncol(d)
	xmat <- matrix(nrow = nr, ncol = nc)

	# Global lowest found value
	globallow <- bestknown
	# Initial solution is empty
	solution <- c()

	# Iteration counters
	branches <- 0
	bounds <- 0		
	ends <- 0

	# Do a heuristic initial guess
	if(verb>=0) print("Performing initial sorting for a good initial guess")
	allowedhclust <- c("complete", "ward", "single", "average", "mcquitty", "median", "centroid")
	if(presort %in% allowedhclust | (presort>=1 & presort<=7)){
		if(is.numeric(presort)) presort <- allowedhclust[presort]
		clust <- hclust(as.dist(d), method=presort)
		d <- d[clust$order, clust$order]	
		returnorder <- order(clust$order)
	}else{
		returnorder <- c(1:nr)	
	}

	if(verb>=0) print("Computing boundaries for minimum distances in possible combinations...")
	# Sorting the whole distance matrix so that we know in what incrementing order the values in a row will be
	# These are required by the 'rowwisesmallestsum'-function
	ranksmat <- t(apply(d, MARGIN=1, FUN=rank, ties.method="first"))
	ordersmat <- t(apply(d, MARGIN=1, FUN=order))
	minsmat <- matrix(nrow=nrow(d), ncol=ncol(d))
	for(i in 1:nrow(d)){
		minsmat[i,] <- unlist(d[i,ordersmat[i,]])
	}
	
	rowwisesmallestsum <- function(
		# Which indices are are free and should be looped through (rows included in minsmat)
		freeindex
	){
		sumsofar <- 0
		for(i in freeindex){
			# Adding g-1 lowest values row-wise (each individual is matched to g-1 other individuals) 
			sumsofar <- sumsofar + sum(minsmat[i,-ranksmat[i,-freeindex]][1:(g-1)], na.rm=TRUE)
		}
		sumsofar		
	}	
		
	progresstemp <- 0
	if(verb>=0) print("Starting branch and bound")
	branch <- function(
		# Binary vector of free indices
		free,
		# Current cost so far
		cost
	){
		if(progresstemp>=progress & verb>=0){
			print(paste("Branching operation #",branches,"..."))
			print(paste("Current best solution: c(", paste(solution[returnorder], collapse=","), ")"))
			print(paste("Current solution cost", globallow))
			print(Sys.time())
			print("")
			progresstemp <<- 0
		}
		progresstemp <<- progresstemp + 1
		branches <<- branches +1
		if(all(!free==0)){
			# At an end node
			if(cost<globallow){
				globallow <<- cost
				solution <<- free
			}
			ends <<- ends + 1
		}else{
			frees <- which(free==0)
			firstindex <- frees[1]
			currentmatch = max(free)+1
			# If only one possibility
			if(currentmatch == nr / g){
				combs <- as.matrix(frees)
			# Else iterate through possibilities
			}else{
				free[firstindex] = currentmatch
				combs = rbind(firstindex, combn(frees[-1], (g-1)))
			}
			for(i in 1:ncol(combs)){
				# Temporary vector of free indexes, reserving some
				freetemp <- free
				freetemp[combs[,i]] = currentmatch
				freestemp = which(freetemp==0)
				# Temporary cost for this new branch
				costtemp <- cost + sum(d[combs[,i], combs[,i]], na.rm=TRUE)
				if(length(freestemp)>0){
					costmin <- costtemp + rowwisesmallestsum(freeindex = freestemp)
				}else{
					costmin <- costtemp
				}
				if(costmin < globallow & branches < maxbranches){
					# Theoretic minimum is better than current global low, we shall explore the node further
					branch(free = freetemp, cost = costtemp)
				}else{
					# Branch is bound
					bounds <<- bounds + 1
				}
			}
		}
	}
	
	# Initiating branch and bound
	branch(
		# Initially no individuals are matched
		free = rep(0, times=nr),
		# Initially the cost is zero
		cost = 0
	)

	solution <- solution[returnorder]

	if(verb>=0){
		print(paste("Branches:", branches))
		print(paste("Bounds:",bounds))
		print(paste("Ends visited:", ends))
		print(paste("Solution cost", globallow))
		print(paste("Solution:",paste(solution,collapse=",")))
	}

	xmat = matrix(0, nrow=nr, ncol=nc)
	uniqsol <- unique(solution)
	for(i in uniqsol){
		indices <- which(solution==uniqsol[i])
		combs <- combn(indices, 2)
		for(j in 1:ncol(combs)){
			xmat[combs[1,j], combs[2,j]] <- xmat[combs[2,j], combs[1,j]] <- 1
		}
	}
	rownames(xmat) <- rownames(d)[returnorder]
	colnames(xmat) <- colnames(d)[returnorder]
	names(solution) <- colnames(d)[returnorder]

	list(branches = branches, bounds = bounds, ends = ends, matrix = xmat, solution = solution, cost = globallow)
}

# Transform a matching vector of form m = {m1, m2, m3, m4, ..., mn} where mn indicate sub-match indices to a binary matching matrix of size n x n
match.vec2mat = function(
	# Vector indicating matched pairs/triplets/...
	# Each unique element is a sub-match
	x
){
	xmat = matrix(0, nrow=length(x), ncol=length(x))
	uniqsol <- unique(x)
	for(i in 1:length(uniqsol)){
		indices <- which(x==uniqsol[i])
		combs <- combn(indices, 2)
		for(j in 1:ncol(combs)){
			xmat[combs[1,j], combs[2,j]] <- xmat[combs[2,j], combs[1,j]] <- 1
		}
	}
	rownames(xmat) <- rownames(x)
	colnames(xmat) <- colnames(x)

	xmat
	# Returning the constructed matrix
}

# Transform a binary matching matrix to a matching vector of form m = {m1, m2, m3, m4, ..., mn} where mn indicate sub-match indices
match.mat2vec = function(
	# Binary matching matrix
	xmat
){
	vec <- rep(0, times=nrow(xmat))
	for(i in 1:length(vec)){
		if(vec[i]==0){
			vec[i] <- max(vec)+1
			vec[which(xmat[i,]==1)] <- vec[i]
		}
	}
	vec
}

# Add averaged dummy observations to the input data matrix in order to create feasible pairs/triplets/quadruplets/...
# Alternatively, add 0-distance 'sinks' to the distance matrix 'd'
# Notice that these two approaches most likely result in different matching solutions
match.dummy = function(
	# Input data matrix
	# where rows equal to individuals and columns to separate biomarkers
	dat,
	# Input (square) distance/dissimilarity matrix
	d,
	# The number of elements in each sub-match, default value is for paired matching
	g = 2
){
	if(!missing(dat)){
		cm = colMeans(dat)
		dummy <- 1
		while(dim(dat)[1]%%g>0)
		{
			dat = rbind(dat, cm)
			rownames(dat)[nrow(dat)] <- paste("Dummy", dummy, sep="")
			dummy <- dummy + 1
		}
		dat
	}else if(!missing(d)){
		if(!dim(d)[1]==dim(d)[2]) stop("Distance matrix d should be square matrix n x n")
		dummy <- 1
		while(!dim(d)[1]%%g==0){
			new.rows <- matrix(0, nrow=1, ncol=ncol(d))
			rownames(new.rows) <- paste("Dummy", dummy, sep="")
			new.cols <- matrix(0, nrow=nrow(d)+1, ncol=1)
			colnames(new.cols) <- paste("Dummy", dummy, sep="")
			d <- cbind(rbind(d, new.rows), new.cols)
			dummy <- dummy + 1
		}
		d
	}else{
		stop("You must provide either data matrix 'dat' or distance/dissimilarity matrix 'd'")
	}
}

# Randomly allocate each member of a sub-match to different groups
match.allocate <- function(
	# Binary matching matrix for the units to be allocated
	# may also be matching vector, which is cast to the matching matrix before allocation
	xmat
){
	# If not matrix, try to cast matching vector to matching matrix
	if(!class(xmat) %in% c("matrix", "data.frame")){
		xmat <- match.vec2mat(xmat)
	}

	# Check that input is in correct format
	if(!dim(xmat)[1] == dim(xmat)[2]) stop("Input should be a square binary matching matrix")
	sums <- apply(xmat, FUN=sum, MARGIN=1)
	if(!all(sums==sums[1])){
		warning("Inconsistent row sums - possibly broken matching matrix")
	}
	sums <- apply(xmat, FUN=sum, MARGIN=2)
	if(!all(sums==sums[1])){
		warning("Inconsistent column sums - possibly broken matching matrix")
	}

	# Perform allocation
	groups <- paste("Group_", LETTERS[1:(sums[1]+1)], sep="")
	groupvec <- vector(length=dim(xmat)[1])
	names(groupvec) <- rownames(xmat)
	diag(xmat) <- 1
	index <- 1
	while(index<dim(xmat)[1]){
		indices <- which(xmat[index,]>0)
		if(length(indices)>0){
			xmat[indices,] <- xmat[,indices] <- 0
			groupvec[indices] <- sample(groups)
		}
		index <- index + 1
	}
	groupvec
}


