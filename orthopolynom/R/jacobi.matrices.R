jacobi.matrices <- function( r )
{
###
### This function returns a list of n real symmetric matrices
### which are the principal minors of the n by n Jacobi matrix
### derived from the the monic recurrence parameters, a and b.
###
### Parameter
### r = a data frame containing the parameters a and b
###
    a <- r$a
    sqrt.b <- sqrt( r$b )
    np1 <- nrow( r )
    n <- np1 - 1
    A <- diag( a[1:n] )
    for ( i in 2:n ) {
        A[i,i-1] <- sqrt.b[i]
        A[i-1,i] <- sqrt.b[i]
    }
    matrices <- as.list( rep( NULL, n ) )
    for ( j in 1:n ) {
        matrices[[j]] <- A[1:j,1:j]
    }
    return( matrices )
}   
