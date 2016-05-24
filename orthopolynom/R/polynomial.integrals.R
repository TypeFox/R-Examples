polynomial.integrals <- function( polynomials )
{
###
###   This function returns a list with n+1 elements
###   containing the integral of the order k polynomials
###   for orders k=0,1,...,n from zero to x
###
###   Parameter
###   polynomials = a list of polynomial objects
###
   require( polynom )
   n <- length( polynomials )
   integrals <- as.list( rep( NULL, n ) )
   j <- 1
   while ( j <= n ) {
      integrals[[j]] <- integral( polynomials[[j]] )
      j <- j + 1
   }
   return( integrals )
}
