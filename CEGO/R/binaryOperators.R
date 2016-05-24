#   Copyright (c) 2014-2015 by Martin Zaefferer, Cologne University of Applied Sciences

###################################################################################
#' Binary String Generator Function
#' 
#' Returns a function that generates random bit-strings of length N.
#' Can be used to create individuals of NK-Landscapes or other problems with binary representation.
#'
#' @param N length of the bit-strings
#'
#' @return returns a function, without any arguments
#'
#' @export
###################################################################################
solutionFunctionGeneratorBinary <- function(N){
	N #lazy evaluation fix, faster than force()
	function()sample(c(0,1),N,replace=TRUE)
}


###################################################################################
#' Bit-flip Mutation for Bit-strings
#' 
#' Given a population of bit-strings, this function mutates all 
#' individuals by randomly inverting one or more bits in each individual. 
#'
#' @param population List of bit-strings
#' @param parameters list of parameters: parameters$mutationRate => mutation rate, specifying number of bits flipped. Should be in range between zero and one
#'
#' @return mutated population
#'
#' @export
###################################################################################
mutationBinary <- function(population, parameters){
  mutationRate <- parameters$mutationRate
	N=length(population[[1]])
	popsize= length(population)
	newpop <- list()
	cmp <- max(min(round(mutationRate*N),N),1)
	for(i in 1:popsize){		
		index<-sample(N,cmp,FALSE,NULL)
		individual <- population[[i]]
		individual[index]=as.numeric(!individual[index])
		newpop <- c(newpop, list(individual))
	}	
	newpop
}

###################################################################################
#' Bit-flip Mutation for Bit-strings (Fast)
#' 
#' Given a population of bit-strings, this function mutates all 
#' individuals by randomly inverting one bit in each individual. 
#' Due to the fixed mutation rate, this is computationally faster.
#'
#' @param population List of bit-strings
#' @param parameters not used
#'
#' @return mutated population
#'
#' @export
###################################################################################
mutationBinaryFast <- function(population, parameters){
	N=length(population[[1]])
	popsize= length(population)
	newpop <- list()
	index <- sample(N,popsize,TRUE,NULL)
	for(i in 1:popsize){		
		individual <- population[[i]]
		individual[index[i]]=as.numeric(!individual[index[i]])
		newpop <- c(newpop, list(individual))
	}	
	newpop
}


###################################################################################
#' Uniform Crossover for Bit Strings
#' 
#' Given a population of bit-strings, this function recombines each
#' individual with another individual by randomly picking bits from each parent. 
#' Note, that \code{\link{optimEA}} will not pass the whole population
#' to recombination functions, but only the chosen parents.
#'
#' @param population List of bit-strings
#' @param parameters not used
#'
#' @return population of recombined offspring
#'
#' @export
###################################################################################
recombinationBinaryUniformCrossoverFast <- function(population, parameters){
	N=length(population[[1]])
	popsize= length(population)
	newpop <- list()	
	for(i in 1:popsize){
		index<-sample(N,N*0.5,FALSE,NULL)
		j <- (1:popsize)[-i][sample(popsize-1,1,FALSE,NULL)] #draw second parent
		parent1 <- population[[i]]
		parent1[-index]=population[[j]][-index] #contribution of second parent			
		newpop <- c(newpop, list(parent1))
	}	
	newpop
}
