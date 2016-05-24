ev.llgp <- function( tab, k )
{
###
### This function returns the subscript of the entering variable at the k-th level
###
### Parameters
### tab = a list of named components that specify the the modified simplex tableau
### k = an integer priority level
###
    sp <- 0
    vsp <- 0.0
    for ( s in 1:tab$nonbasics ) {
        if ( tab$ti[k,s] > 0.0 ) {
            if ( ( k == 1 ) || ( neg.ind.rows(tab, k, s) == 0 ) ) {
                if ( tab$ti[k,s] > vsp ) {
                    sp <- s
                    vsp <- tab$ti[k,s]
                }
            }
        }
    }
    return( sp )
}
