neg.ind.rows <- function( tab, k, s )
{
###
### This function returns a count of the number of negative index values above I[k,s]
###
### Parameters
### tab = a list of named components that specify the modified simplex tableau
### k = an integer priority level
### s = an integer index for non-basic variable
###
    if ( k == 1 )
        return( 0 )
    n <- 0
    l <- 1
    while ( l < k ) {
        if ( tab$ti[l,s] < 0.0 )
            n <- n + 1
        l <- l + 1
    }
    return( n )
}
