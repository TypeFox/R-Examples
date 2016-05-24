#   Copyright (c) 2014-2015 by Martin Zaefferer, Cologne University of Applied Sciences

###################################################################################
#' NK-Landscape Benchmark Creation
#'
#' Function that generates a NK-Landscapes.
#'
#' @param N length of the bit strings
#' @param K number of neighbours contributing to fitness of one position
#' @param PI vector, giving relative positions of each neighbour in the bit-string
#' @param g set of fitness functions for each possible combination of string components. Will be randomly determined if not specified. Should have N rows, and 2^(K+1) columns.
#'
#' @return the function of type cost=f(bitstring). Returned fitness values will be negative, for purpose of minimization.
#'
#' @examples
#' fun <- benchmarkGeneratorNKL(6,2)
#' fun(c(1,0,1,1,0,0))
#' fun(c(1,0,1,1,0,1))
#' fun(c(0,1,0,0,1,1))
#' fun <- benchmarkGeneratorNKL(6,3)
#' fun(c(1,0,1,1,0,0))
#' fun <- benchmarkGeneratorNKL(6,2,c(-1,1))
#' fun(c(1,0,1,1,0,0))
#' fun <- benchmarkGeneratorNKL(6,2,c(-1,1),g=matrix(runif(48),6))
#' fun(c(1,0,1,1,0,0))
#' fun(sample(c(0,1),6,TRUE))
#' 
#' @export
###################################################################################
benchmarkGeneratorNKL <- function(N=10,K=1,PI=1:K,g=NULL){  
	if(is.null(g)){ # generate the fitness subfunctions
		g <- matrix(runif(N*2^(K+1)),N)
	}
	bits = 2^(0:K)
	bits
	g
	N
	K
	PI
	function(x){
		usum=0
		for(i in 1:N){
			xx <- x[c(i,((i+PI-1)%%N)+1)] #select current and impacting (neighbouring) bits (circular)
			usum=usum+ g[i,sum(bits*xx)+1]
		} 
		-usum/N #minus for minimization
	}
}