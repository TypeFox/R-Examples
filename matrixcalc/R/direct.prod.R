direct.prod <- function( x, y )
{
###
### This function computes the direct product of two vectors or
### matrices resulting in a matrix
###
### Arguments
### x = a numeric matrix or vector
### y = a numeric matrix or vector
###
    if ( is.vector( x ) ) {
        A <- matrix( x )
    }
    else {
        if ( is.matrix( x ) ) {
            A <- x
        }
        else {
            stop( "Argument x is not a matrix or vector" )
        }
    }
    if ( is.vector( y ) ) {
        B <- matrix( y )
    }
    else {
        if ( is.matrix( y ) ) {
            B <- y
        }
        else {
            stop( "Argument y is not a matrix or vector" )
        }
    }
###
### construct a list of lists.  The number of elements in the higher level
### is the number of rows in A.  The number of elements in the lower level
### is the number of columns in A.  The value for higher level index i
### and lower level index j is A[i,j] * B
###
    matrices <- list()
    for ( i in 1:nrow( A ) ) {
        matrices[[i]] <- list()
        for ( j in 1:ncol( A ) ) {
            matrices[[i]][[j]] <- A[i,j] * B
        }
    }
###
### construct a list of matrices.  The number of elements is the number
### of rows in A.  The value for element i is the column binding of all the elements
### in the higher level matrices[[i]].
###
    row.matrices <- list()
    for ( i in 1:nrow( A ) ) {
        row.matrices[[i]] <- matrices[[i]][[1]]
        if ( ncol( A ) > 1 ) {
            for ( j in 2:ncol( A ) ) {
                row.matrices[[i]] <- cbind( row.matrices[[i]], matrices[[i]][[j]] )
            }
        }
    }
###
### construct the result matrix C as the row binding of all of the elements in row.matrices
###
    C <- row.matrices[[1]]
    if ( length( row.matrices ) > 1 ) {
        for ( i in 2:length( row.matrices ) ) {
            C <- rbind( C, row.matrices[[i]] )
        }
    }    
    return( C )
}
