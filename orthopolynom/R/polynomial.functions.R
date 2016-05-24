polynomial.functions <- function( polynomials, ... )
{
###
### this function returns a list of functions one for each of the
### polynomials in the arguments
###
### Parameters
### polynomials = a list of polynomial objects
### ... = further arguments to be passed to or from methods
###
    functions <- list()
    n <- length( polynomials )
    for ( j in 1:n ) {
        functions[[j]] <- as.function( polynomials[[j]], ... )
    }
    return( functions )
}
