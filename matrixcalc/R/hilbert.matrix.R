hilbert.matrix <- function(n) 
{   
###
### this function returns an n by n Hilbert matrix
###
### Parameter
### n = the row (column) dimension of the matrix
###
    if ( n != trunc( n ) )
        stop( "argument n is not an integer" )
    if ( n < 2 )
        stop( "argument n is less than 2" )
    i <- 1:n
    X <- 1 / outer(i - 1, i, "+")
    return( X  )
}
