symmetric.pascal.matrix <- function( n )
{
###
### this function returns an n by n Pascal matrix
###
### Parameter
### n = the row( column ) dimension of the matrix
###
    if ( n <= 0 )
        stop( "argument n is not positive" )
    if ( n != trunc( n ) )
        stop( "argument n is not an integer" )
    nm1 = n-1
    n.over.r <- function(n, r) { prod(1:n) / (prod(1:(n-r)) * prod(1:r) ) }
    X <- rep(1, nm1)
    for ( i in 1:nm1 )
        for ( j in 1:nm1 )
            X <- c(X, n.over.r(i+j, j))
    X <- cbind(rep(1, nm1+1), matrix(X, byrow = TRUE, ncol = nm1))
    return( X  )
}
