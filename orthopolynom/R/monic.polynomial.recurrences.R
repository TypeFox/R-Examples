monic.polynomial.recurrences <- function( recurrences )
{
###
### This function returns a data frame with parameters
### required to construct monic orthogonal polynomials
### based on the standard recurrence relation for
### the non-monic polynomials
###
### Parameter
### recurrences = the data frame of recurrence parameters c, d, e and f
###
    np1 <- nrow( recurrences )
    n <- np1 - 1
    monic.recurrences <- data.frame( matrix( nrow=np1, ncol=2 ) )
    names( monic.recurrences ) <- c( "a", "b" )
    c <- recurrences$c
    d <- recurrences$d
    e <- recurrences$e
    f <- recurrences$f
    j <- 0
    k <- 1
    while( j <= n ) {
        if ( e[k] == 0 )
            monic.recurrences[k,"a"] <- 0
        else
            monic.recurrences[k,"a"] <- - ( d[k] / e[k] )
        if ( j == 0 )
            monic.recurrences[k,"b"] <- 0
        else {
            if ( e[k-1] == 0 || e[k] == 0 )
                monic.recurrences[k,"b"] <- 0
            else
                monic.recurrences[k,"b"] <- ( c[k-1] * f[k] ) / ( e[k-1] * e[k] )
        }
        j <- j + 1
        k <- k + 1
    }
    return( monic.recurrences )
}
