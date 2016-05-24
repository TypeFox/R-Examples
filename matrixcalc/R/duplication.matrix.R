duplication.matrix <- function( n=1 )
{
###
### this function returns a matrix sith n * n rows and n * (n + 1 ) / 2
### columns that transforms vech( A ) to vec( A ) where A is a symmetric n by n matrix
###
### Parameter
### n = the order of the matrix
###
    return( D.matrix( n ) )
}
