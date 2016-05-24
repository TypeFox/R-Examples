frobenius.prod <- function( x, y )
{
###
### this function calculates the Frobenius inner product product of two matrices x and y.
### the matrices must the same row and column order
###
### Parameters
### x = a numeric matrix or vector object
### y = a numeric matrix or vector object
###
    return( sum( hadamard.prod(x, y) ) )
}
