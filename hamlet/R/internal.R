#
# Internal functions for use in the hamlet-package
# At the time of writing this they are mainly different functions used in the Genetic Algorithm (GA)
#


.ga.breed <- function(x, y, d){
	# Find out a row in d*x (distance mat * matching mat) where difference is greatest in parents and adjust towards better
	rowdifs <- unlist(apply(x*d, MARGIN=1, FUN=sum)) - unlist(apply(y*d, MARGIN=1, FUN=sum))
	ord <- which.min(rowdifs)
	w <- match.vec2mat(x)[ord,] - match.vec2mat(y)[ord,]
	ind <- 1:length(w)
	w1 <- which(w>0) #w1 <- which(w==1)
	w2 <- which(w<0) #w2 <- which(w==-1)
	ind[w1] <- w2
	ind[w2] <- w1
	y[ind]
}

.ga.mutate = function(x){
	# A point mutation in the solution vector is done by randomly permutating a row/column in matching matrix or two elements in the solution vector
	ran <- sample(1:length(x), 2, replace=F)
	temp <- x[ran[1]]; x[ran[1]] <- x[ran[2]]; x[ran[2]] <- temp;
	x
}	
.ga.fitness = function(x, d){
		uniqs <- unique(x) # All unique submatches in the solution
		# Take all distances from 'd' according to what is matched and sum to obtain the target function value
		# Sum over all submatches
		sum(unlist(lapply(uniqs, FUN=function(z){
			w <- which(x == z)
			# Sum over pairs of individuals within a submatch
			sum(unlist(apply(combn(w, 2), MARGIN=2, FUN=function(y) d[y[1], y[2]])))
		}))) *2 # Multiply by two due to symmetricity in the distance/matching matrices
	}	
.ga.weight = function(fitnesses){
		#w <- rank(fitnesses)^10
		w <- rank(fitnesses)^100
		w/sum(w)
	}
.ga.step = function(pops, fitnesses, weights, nmutate, ndeath, mutate, breed, fitness, weight, d){
		if(missing(weights)) weights <- weight(fitnesses)
		# Mutations are sampled according to fitness weights with replacement
		indexmutate <- sample(1:ncol(pops), prob = weights, size = nmutate, replace = T)
		# Solution deaths are sampled without replacement according to fitness weights
		indexdeath <- sample(1:ncol(pops), prob = weights, size = ndeath, replace = F)
		# The dual parents of a new solution to replace a dead one are sampled according to an inverse fitness weight
		indexparents <- replicate(n = ndeath, expr = sample(1:ncol(pops), size=2, prob=1-weights, replace=F))
		
		# Create new generation of solutions
		# Mutations
		# bad solutions (relative to the rest) are more likely to receive mutations
		for(i in 1:length(indexmutate)){
			pops[,indexmutate[i]] <- mutate(pops[,indexmutate[i]])
		}
		# Every dying solution is replaced by an offspring produced by two existing solutions
		# better solutions are more likely to be parents of an offspring
		# Create offspring
		offspring <- lapply(1:length(indexdeath), FUN=function(z){
			breed(pops[,indexparents[1,z]], pops[,indexparents[2,z]], d = d)
		})
		# Place offspring to the population matrix (do in two steps so we don't replace previous entries before breeding is done)
		for(i in 1:length(offspring)){
			pops[,indexdeath[i]] <- offspring[[i]]
		}
		# Return new generation
		fitnesses[union(indexmutate, indexdeath)] <- unlist(apply(pops[,union(indexmutate, indexdeath)], MARGIN=2, FUN=function(z) fitness(z, d = d)))
		weights <- weight(fitnesses)
		list(pops = pops, fitnesses = fitnesses, weights = weights)
	}	
.ga.init = function(popsize = popsize, g, d){
		replicate(n = popsize, expr = sample(rep(1:(nrow(d)/g), each=g)))
	}



