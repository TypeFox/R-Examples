dv.llgp <- function( tab, sp )
{
###
### This function returns the index of the departing variable from the basis
### for a lexicographical linear goal programming problem
###
### Parameters
### tab = a list of named components that specify the modified simplex tableau
### sp = integer index of the non-basic variable
###
    ip <- 0
    for ( i in 1:tab$objectives ) {
        if ( tab$te[i,sp] > 0.0 ) {
            v <- tab$tb[i] / tab$te[i,sp]
            if ( ip == 0 ) {
                vip <- v
                ip <- i
            }
            else if ( v < vip ) {
                vip <- v
                ip <- i
            }
            else if ( v == vip ) {
                ip <- dv.tie( tab, i, ip )
            }
        }    
    }
    return( ip )
}
