# Copyright (C) 2011 Jelmer Ypma. All Rights Reserved.
# This code is published under the GPL.
#
# File:   SparseGridGetSeq.R
# Author: Jelmer Ypma
# Date:   31 March 2011
#
# This code is based on code from www.sparse-grids.de
# with permission from the authors.
#
# Description: function for generating matrix of all rowvectors in N^D_{norm}.
#			   for internal use.
#
# Input: 
#		dimension = dimension, will be number of columns in output
#		norm	  = row sum of elements each of the rows has to have
#
# Output: 
#		a matrix with dimension columns. Each row represents one vector with
#		all elements >= 1 and the sum of elements == norm
#		

SparseGridGetSeq <- function( dimension, norm ) {
    seq.vec       <- rep(0, dimension)
    a             <- norm - dimension
    seq.vec[ 1 ]  <- a
    fs            <- matrix( seq.vec, nrow=1, ncol=length(seq.vec) )
    cnt           <- 1
    while ( seq.vec[ dimension ] < a ) {
        if ( cnt == dimension ) {
            for (i in seq(cnt-1, 1, -1) ) {
                cnt <- i
                if ( seq.vec[ i ] != 0 ) { break }
            }
        }
        seq.vec[ cnt ]  <- seq.vec[ cnt ] - 1
        cnt             <- cnt + 1
        seq.vec[ cnt ]  <- a - sum( seq.vec[1:(cnt-1)] )
        if ( cnt < dimension ) {
            seq.vec[ (cnt+1):dimension ] <- rep(0, dimension - cnt)
        }
        fs = rbind( fs, seq.vec )
    }
    fs <- fs + 1
    
    return( fs )
}
