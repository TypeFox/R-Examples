polynomial.derivatives <- function( polynomials )
{
###
###   This function returns a list with n+1 elements
###   containing the derivatives of the order k polynomials
###   for orders k=0,1,...,n
###
###   Parameter
###   polynomials = a list of polynomial objects
###
   require( polynom )
   n <- length( polynomials )
   derivatives <- as.list( rep( NULL, n ) )
   j <- 1
   while ( j <= n ) {
      derivatives[[j]] <- deriv( polynomials[[j]] )
      j <- j + 1
   }
   return( derivatives )
}
