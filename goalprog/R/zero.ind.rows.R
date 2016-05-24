zero.ind.rows <- function( tab, s )
{
###
### This function returns the number of zero values in column s
###
### Parameters
### tab = a list of named components that specifies the modified simplex tableau
### s = an integer index for a non-basic variable
###
    n <- 0
    for ( k in 1:tab$levels ) {
        if ( tab$ti[k,s] == 0.0 )
            n <- n + 1
    }
    return( n )
}
