pos.ind.rows <- function( tab, k, s )
{
###
### this function returns the number of positive index values above I(k,s)
###
### Parameters
### tab = a list of named components with the modified simplex tableau
### k = an integer priority level
### s = an integer index for a non-basic variable
###
    if ( k == 1 )
        return( 0 )
    n <- 0
    l <- 1
    while( l < k ) {
        if (tab$ti[l,s] > 0.0 )
            n <- n + 1
        l <- l + 1    
    }
    return( n )
}
