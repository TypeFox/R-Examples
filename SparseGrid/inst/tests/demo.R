# Copyright (C) 2011 Jelmer Ypma. All Rights Reserved.
# This code is published under the GPL.
#
# File:   demo.R
# Author: Jelmer Ypma
# Date:   31 March 2011
#
# This code is based on the example from www.sparse-grids.de
# with permission from the authors.

library('SparseGrid')

# set seed of random number generator
set.seed( 3141 )

dimension   <- 10   # dimensions
maxk        <- 4    # max. accuracy level (pol. exactness wil be 2k-1)

# integrand: some function that evaluates g(x): (R times D)->(R times 1)
func <- function( x, sigma=2 ) {
    return( apply( exp(-.5*(x/sigma)^2)/sqrt(2*pi*sigma^2), 1, prod ) )
} 

# calculate "correct" result of integration between 0 and 1:
trueval <- (pnorm(1, sd=2) - pnorm(0, sd=2))^dimension

# create matrix to hold results
res <- matrix( NA, nrow=maxk-1, ncol=5 )
colnames( res ) <- c("D", "k", "nodes", "SG error", "Sim. error")
rownames( res ) <- rep( "", maxk-1 )

# loop over different accuracy levels
for ( k in 2:maxk ) {

    # sparse grids integration
    tmp.grid <- createSparseGrid('KPU', dimension, k)
    x <- tmp.grid$nodes
    w <- tmp.grid$weights
  	g <- func( x )
    SGappr <- sum(w * g)
    SGerror <- sqrt((SGappr - trueval)^2) / trueval
    
    # simulation with the same number of nodes, 1000 simulation repetitions
    numnodes <- length(w)
    sim <- rep(0, 1000)
    for (r in 1:1000) {
        x <- matrix( runif( numnodes * dimension ), nrow=numnodes, ncol=dimension )
        g <- func( x )
        sim[ r ] <- mean( g )       # is sum(w * g) where weights are 1/numnodes
    }
    Simerror = sqrt(mean((sim-trueval)^2)) / trueval
    
    # save results in row of matrix res
    res[k-1,] <- c(dimension, k, numnodes, SGerror, Simerror)
}
res
