quadrature.rules <- function( recurrences, inner.products )
{
###
### This function returns a list with n elements
### containing the order k quadrature rule data frames
### for orders k=1,2,...n.  An order k quadrature data frame
### contains the roots and abscissa values for the degree k polynomial
###
### Parameters
### recurrences = a data frame with recurrence parameters c, d, e and f
### inner.products = a vector of norm squared values
###
    np1 <- nrow( recurrences )
    n <- np1 - 1
    rules <- as.list( rep( NULL, n ) )
    monic.recurrences <- monic.polynomial.recurrences( recurrences )
    matrices <- jacobi.matrices( monic.recurrences )
    matrix.eigens <- lapply( matrices, eigen )
    roots <- polynomial.roots( monic.recurrences )
    h.0 <- inner.products[1]
    for ( k in 1:n ) {
        values <- matrix.eigens[[k]]$values
        vectors <- matrix.eigens[[k]]$vectors
        x <- values
        w <- rep( 0, k )
        for ( j in 1:k ) {
            v.j <- vectors[1,j]
            w[j] <- h.0 * v.j * v.j
        }
        rule <- data.frame( cbind( x, w ) )
        names( rule ) <- c( "x", "w" )
        rules[[k]] <- rule
    }
    return( rules )
}
