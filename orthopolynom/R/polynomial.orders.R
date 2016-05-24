polynomial.orders <- function( polynomials )
{
#
### This function returns a vector with n elements
### containing the orders of the polynomials in the list
###
### Parameter
### polynomials = a list of polynomial objects
###
    require( polynom )
    n <- length( polynomials )
    orders <- rep( 0, n ) 
    j <- 1
    while ( j <= n ) {
        coefficients <- as.vector( polynomials[[j]] )
        orders[j] <- length( coefficients ) - 1
        j <- j + 1
    }
    return ( orders )
}   
