polynomial.powers <- function( polynomials )
{
###
### This function returns a list of vectors that define linear combinations
### of the given polynomials equal to powers of the polynomial argument.
###
### Parameters:
### polynomials = a list of polynomial objects in increase polynomial order
###
    require( polynom )
    coefficients <- polynomial.coefficients( polynomials )
    np1 <- length( polynomials )
    powers <- as.list( rep( NULL, np1 ) )
    lambda.matrix <- matrix( 0, nrow=np1, ncol=np1 )
    for ( j in 1:np1 ) {
        lambda <- coefficients[[j]]
        for ( k in 1:j ) {
            lambda.matrix[j,k] <- lambda[k]
        }
    }
    alpha.matrix <- solve( lambda.matrix )
    for ( j in 1:np1 ) {
        alpha <- rep( 0, j )
        for ( k in 1:j ) {
            alpha[k] <- alpha.matrix[j,k]
        }
        powers[[j]] <- alpha
    }
    return( powers )
}
