# Copyright (C) 2011 Jelmer Ypma. All Rights Reserved.
# This code is published under the GPL.
#
# File:   createMonteCarloGrid.R
# Author: Jelmer Ypma
# Date:   19 November 2011
#
# Description: function to calculate nodes and weights for numerical integration 
#              using Monte Carlo simulations.
#
# Input: 
#     rng       : function that generates random numbers. The first argument of this function
#                 should be called 'n'. Examples are the R built-in functions rnorm and runif 
#                 for random numbers from a standard normal or uniform distribution.
#     dimension : dimension of the integration problem.
#     num.sim   : number of simulated integration nodes.
#     ...       : arguments that will be passed to the random number generator rng.
#
# Output: 
#     list with two elements:
#		nodes    	= matrix of nodes with dim columns 
#     	weights    = row vector of corresponding weights
#
# Examples:
#	set.seed( 3141 )
#	createMonteCarloGrid( runif, dimension=2, num.sim=10 )
#	createMonteCarloGrid( rnorm, dimension=2, num.sim=10 )
#	createMonteCarloGrid( rnorm, dimension=2, num.sim=10, mean=2, sd=5 )

createMonteCarloGrid <- function( rng, dimension, num.sim, ... ) {
	if ( !is.function( rng ) ) {
		stop('argument rng in createMonteCarloGrid should be a function.')
	}
	
	if ( !names( formals( rng ) )[1] == "n" ) {
		stop('first argument of supplied rng is not called "n". Examples of random number generators that you can use are rnorm and runif.')
	}
	
	return(
		list( "nodes" = matrix( rng( num.sim * dimension, ... ), nrow=num.sim, ncol=dimension ),
		      "weights" = rep( 1 / num.sim, num.sim ) 
		)
	)
}
