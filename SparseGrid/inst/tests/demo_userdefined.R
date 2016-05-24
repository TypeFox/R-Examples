# Copyright (C) 2011 Jelmer Ypma. All Rights Reserved.
# This code is published under the GPL.
#
# File:   demo_userdefined.R
# Author: Jelmer Ypma
# Date:   31 March 2011
#
# This example shows how to use a user-defined function
# to create a sparse grid, and shows what can go wrong

library('SparseGrid')   # USES createSparseGrid
require('statmod')      # USES guass.quad

ghq <- function(n) 
{ 
    grid <- gauss.quad( n, 'hermite' )
    return( 
        list(
            "nodes" = grid$nodes * sqrt(2), 
            "weights" = grid$weights / sqrt(pi) 
            )
          )
}

ghq.symmetric <- function(n)
{
	grid <- gauss.quad( n, 'hermite' )
	select.positive <- grid$nodes >= 0
	return( 
        list(
            "nodes" = grid$nodes[ select.positive ] * sqrt(2), 
            "weights" = grid$weights[ select.positive ] / sqrt(pi) 
            )
          )
}

#
res1a <- createSparseGrid(ghq, 1, 5)
res1b <- createSparseGrid(ghq.symmetric, 1, 5, sym=TRUE)	# misses one node, which is wrong
res1c <- createSparseGrid('GQN', 1, 5)

# there is a duplicate node in res2a, because of rounding errors
# 
res2a <- createSparseGrid(ghq, 2, 3)						# has one duplicate node, which is slightly inefficient
res2b <- createSparseGrid(ghq.symmetric, 2, 3, sym=TRUE)	# misses one node, which is wrong
res2c <- createSparseGrid('GQN', 2, 3)
