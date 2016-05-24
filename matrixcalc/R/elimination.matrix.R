elimination.matrix <- function( n )
{
###
### this function returns a matrix sith n * n comuns and n * (n + 1 ) / 2
### rows that transforms vec( A ) to vech( A ) where A is a symmetric n by n matrix
###
### Parameter
### n = the order of the matrix
###
    return( L.matrix( n ) )
}
