#   Copyright (c) 2014-2015 by Martin Zaefferer, Cologne University of Applied Sciences

#TODO:
# insert/shift mutation
# reversal mutation

###################################################################################
#' Permutation Generator Function
#' 
#' Returns a function that generates random permutations of length N.
#' Can be used to generate individual solutions for permutation problems, e.g., Travelling Salesperson Problem
#'
#' @param N length of the permutations returned
#'
#' @return returns a function, without any arguments
#'
#' @export
###################################################################################
solutionFunctionGeneratorPermutation <- function(N){
	N #lazy evaluation fix, faster than force()
	function()sample(1:N,replace=FALSE)
}

###################################################################################
#' Interchange Mutation for Permutations
#' 
#' Given a population of permutations, this function mutates all 
#' individuals by randomly interchanging two arbitrary elements of the permutation.
#'
#' @param population List of permutations
#' @param parameters list of parameters, currently only uses parameters$mutationRate, 
#' which should be between 0 and 1 (but can be larger than 1). The mutation rate determines the number of interchanges
#' performed, relative to the permutation length (N). 0 means none. 1 means N interchanges.
#' The default is 1/N.
#'
#' @return mutated population
#'
#' @export
###################################################################################
mutationPermutationInterchange <- function(population, parameters=list()){
  N <- length(population[[1]])	
  if(is.null(parameters$mutationRate)) parameters$mutationRate <- 1/N
	mrate <- parameters$mutationRate
	popsize <- length(population)
	mutations <- ceiling(N * mrate)
	if(mutations==0)
		return(population)
	samples <- mutations * popsize
	newpop <- list()
	index1 <- sample.int(N,samples,TRUE,NULL) 
	index2 <- sample.int(N,samples,TRUE,NULL)
	for(i in 1:popsize){				
		individual <- population[[i]]
		if(mutations == 1){
			val1= individual[index1[i]]
			individual[index1[i]]= individual[index2[i]]
			individual[index2[i]]= val1
		}else{
			j <- ((i-1)*mutations+1) : (i*mutations)
			for(jj in j){
				i1 <- index1[jj]
				i2 <- index2[jj]
				val1= individual[i1]
				individual[i1]= individual[i2]
				individual[i2]= val1
			}			
		}
		newpop <- c(newpop, list(individual))
	}	
	newpop
}
#mutationPermutationInterchange(list(1:5),list(mutationRate=1))

###################################################################################
#' Swap Mutation for Permutations
#' 
#' Given a population of permutations, this function mutates all 
#' individuals by randomly interchanging two adjacent elements of the permutation.
#'
#' @param population List of permutations
#' @param parameters list of parameters, currently only uses parameters$mutationRate, 
#' which should be between 0 and 1 (but can be larger than 1). The mutation rate determines the number of swaps
#' performed, relative to the permutation length (N). 0 means none. 1 means N swaps.
#' The default is 1/N.
#'
#' @return mutated population
#'
#' @export
###################################################################################
mutationPermutationSwap <- function(population,parameters=list()){
	N <- length(population[[1]])
  if(is.null(parameters$mutationRate)) parameters$mutationRate <- 1/N
	mrate <- parameters$mutationRate
	popsize <- length(population)
	mutations <- ceiling(N * mrate)
	if(mutations==0)
		return(population)
	samples <- mutations * popsize
	newpop <- list()
	index1 <- sample.int(N-1,samples,TRUE,NULL) 
	index2 <- index1 +1
	for(i in 1:popsize){				 # TODO: after index2, the same code as in interchange mutation is used. code should be merged.
		individual <- population[[i]]
		if(mutations == 1){
			val1 <- individual[index1[i]]
			individual[index1[i]] <- individual[index2[i]]
			individual[index2[i]] <- val1
		}else{
			j <- ((i-1)*mutations+1) : (i*mutations)
			for(jj in j){
				i1 <- index1[jj]
				i2 <- index2[jj]
				val1= individual[i1]
				individual[i1]= individual[i2]
				individual[i2]= val1
			}	
		}		
		newpop <- c(newpop, list(individual))
	}
	newpop
}
#mutationPermutationSwap(list(1:5),list(mutationRate=1))


###################################################################################
#' Reversal Mutation for Permutations
#' 
#' Given a population of permutations, this function mutates all 
#' individuals by randomly selecting two indices, and reversing the respective sub-permutation.
#'
#' @param population List of permutations
#' @param parameters list of parameters, currently only uses parameters$mutationRate, 
#' which should be between 0 and 1 (but can be larger than 1). The mutation rate determines the number of reversals
#' performed, relative to the permutation length (N). 0 means none. 1 means N reversals.
#' The default is 1/N.
#'
#' @return mutated population
#'
#' @export
###################################################################################
mutationPermutationReversal <- function(population, parameters=list()){
	N <- length(population[[1]])
  if(is.null(parameters$mutationRate)) parameters$mutationRate <- 1/N  
	mrate <- parameters$mutationRate
	popsize <- length(population)
	mutations <- ceiling(N * mrate)
	if(mutations==0)
		return(population)
	samples <- mutations * popsize
	newpop <- list()
	index1 <- sample.int(N,samples,TRUE,NULL) 
	index2 <- sample.int(N,samples,TRUE,NULL)
	for(i in 1:popsize){				
		individual <- population[[i]]
		if(mutations == 1){
			individual[index1[i]:index2[i]] <- individual[index2[i]:index1[i]]
		}else{
			j <- ((i-1)*mutations+1) : (i*mutations)
			for(jj in j){
				i1 <- index1[jj]
				i2 <- index2[jj]
				individual[i1:i2] <- individual[i2:i1]
			}	
		}		
		newpop <- c(newpop, list(individual))
	}	
	newpop
}
#mutationPermutationReversal(list(1:5),list(mutationRate=1))


###################################################################################
#' Cycle Crossover (CX) for Permutations
#' 
#' Given a population of permutations, this function recombines each
#' individual with another random individual.
#' Note, that \code{\link{optimEA}} will not pass the whole population
#' to recombination functions, but only the chosen parents.
#'
#' @param population List of permutations
#' @param parameters not used
#'
#' @return population of recombined offspring
#'
#' @export
###################################################################################
recombinationPermutationCycleCrossover <- function(population, parameters){
	popsize <- length(population)
	newpop <- list()	
	for(i in 1:popsize){
		j <- (1:popsize)[-i][sample.int(popsize-1,1,FALSE,NULL)] #draw second parent
		parent1 <- population[[i]]
		parent2 <- population[[j]]
		e1 <- parent1[1]
		e2 <- parent2[1]		
		parent2[1] <- e1
		while(e1 != e2){
			e1 <- e2
			rplc <- which(parent1==e1)
			e2 <- parent2[rplc]
			parent2[rplc] <- e1
		}		
		newpop <- c(newpop, list(parent2))
	}	
	newpop
}