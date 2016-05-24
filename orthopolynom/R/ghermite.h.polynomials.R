ghermite.h.polynomials <- function( n, mu, normalized=FALSE )
{
###
###   This function returns a list with n+1 elements
###   containing the order k generalized Hermite polynomials, Hk(mu,x),
###   for orders k=0,1,...,n
###
###   Parameters
###   n = integer highest polynomial order
###   mu = the polynomial parameter
###   normalized = a boolean value.  if true, the polynomials are normalized
###
   recurrences <- ghermite.h.recurrences( n, mu, normalized )
   if ( normalized ) {
       h.0 <- gamma( mu + 0.5 )
       p.0 <- polynomial( c( 1 / sqrt( h.0 ) ) )
       polynomials <- orthonormal.polynomials( recurrences, p.0 )
       return( polynomials )
   }
   else {
       polynomials <- orthogonal.polynomials( recurrences )
       return( polynomials )
   }    
   return( NULL )
}
