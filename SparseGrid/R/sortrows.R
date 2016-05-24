# Copyright (C) 2011 Jelmer Ypma. All Rights Reserved.
# This code is published under the GPL.
#
# File:   sortrows.R
# Author: Jelmer Ypma
# Date:   31 March 2011
#
# Description:
#			   for internal use.
# Input: 
#
# Output: 
#

# order uses a stable sort algorithm
# http://www.mail-archive.com/r-help@r-project.org/msg57138.html
sortrows <- function( A, index.return=FALSE ) {

    if( !is.matrix( A ) ) {
        stop('SparseGrid:::sortrows expects a matrix as argument A.')
    }
    
    A.nrow <- nrow( A )
    A.ncol <- ncol( A )
    
    if ( index.return==TRUE ) {
        indices <- 1:nrow( A )
    }

    for ( col.cnt in seq( ncol(A), 1, -1 ) ) {
        tmp.indices <- order(A[,col.cnt])

        # drop=FALSE prevent dropping the matrix dimension.
        # We need the return value to be a matrix.
        # If we have only one row, then A[1,] (without drop=FALSE) 
        # returns a vector.
        # Example:
        #       > A <- matrix( c(1,2,3), nrow=1, ncol=3 )
        #       > A
        #            [,1] [,2] [,3]
        #       [1,]    1    2    3
        #       > A[1,]                 # drops the matrix dimensions
        #       [1] 1 2 3
        #       > A[1,,drop=FALSE]      # keeps the dimensions intact
        #            [,1] [,2] [,3]
        #       [1,]    1    2    3
        A <- A[ tmp.indices, , drop=FALSE ]
        
        if ( index.return==TRUE ) {
            indices <- indices[ tmp.indices ]
        }
    }
    
    if ( index.return==TRUE ) {
        res <- list("x" = matrix( A, nrow=A.nrow, ncol=A.ncol ), "ix" = indices )
    }
    else {
        res <- matrix( A, nrow=A.nrow, ncol=A.ncol )
    }
    return( res )
}
