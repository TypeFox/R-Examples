polynomial.values <- function( polynomials, x )
{
###
### This function returns a list with n+1 elements
### containing the values of the order k polynomials
### for orders k=0,1,...,n and for the given argument x
###
### Parameters
### polynomials = a list of polynomial objects
### x = the argument which can be any numeric object
###
    require( polynom )
    n <- length( polynomials )
    values <- as.list( rep( NULL, n ) )
    j <- 1
    while ( j <= n ) {
        values[[j]] <- predict( polynomials[[j]], x )
        j <- j + 1
    }
    return( values )
}
