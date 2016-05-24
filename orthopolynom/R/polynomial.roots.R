polynomial.roots <- function( m.r )
{
###
###   This function returns a list with n elements
###   containing the roots of the order k polynomials
###   for orders k=1,...,n using a data frame with
###   the monic polynomial recurrence parameters a and b
###
###   Parameter
###   m.r = monic recurrence data frame with parameters a and b
###
   matrices <- jacobi.matrices( m.r )
   n <- length( matrices )
   eigen.list <- lapply( matrices, eigen )
   roots <- as.list( rep( NULL, n+1 ) )
   j <- 1
   while ( j <= n ) {
      roots[[j+1]] <- eigen.list[[j]]$values
      j <- j + 1
   }
   return( roots )
}
