polynomial.coefficients <- function( polynomials )
{
#
### This function returns a list with n elements
### containing the vector of coefficients for the polynomials
###
### Parameter
### polynomials = a list of polynomial objects
###
    require( polynom )
    n <- length( polynomials )
    coefficients <- as.list( rep( NULL, n ) )
    j <- 1
    while ( j <= n ) {
        coefficients[[j]] <- as.vector( polynomials[[j]] )
        j <- j + 1
    }
    return ( coefficients )
}   
