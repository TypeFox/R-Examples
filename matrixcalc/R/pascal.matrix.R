pascal.matrix <- function( n )
{
###
### this function returns an n by n Pascal matrix
###
### Parameter
### n = the row( column ) dimension of the matrix
###
    S <- symmetric.pascal.matrix( n )
    luS <- lu.decomposition( S )
    P <- luS$L
    return( P )
}
