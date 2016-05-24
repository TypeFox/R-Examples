dv.tie <- function( tab, i, ip )
{
###
### This function returns an index that resolve ties for departing variables
###
### Parameters
### tab = a list of named components that specify the modified simplex tableau
### i = an integer index for a basic variable
### ip = an integer index for a basic variable
###
    for ( k in 1:tab$levels ) {
        if ( tab$tu[i,k] < tab$tu[ip,k] )
            return ( ip )
        else if ( tab$tu[i,k] > tab$tu[ip,k] ) {
            return ( i )
        }    
    }
    return ( i )
}
