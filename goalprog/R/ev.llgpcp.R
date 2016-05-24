ev.llgpcp <- function( tab, k )
{
###
### This function returns the subscript of the entering variable at the k-th level
### taking into account complementary pivoting.
###
### Parameters
### tab = a list of named components that specify the augmented modified simplex tableau
### k = an integer priority level
###
    sp <- 0
    vsp <- 0.0
###
### determine if complementary variable is in the bais
###
    for ( s in 1:tab$nonbasics ) {
        if ( check.ev.cp( tab, s ) ) {
            if ( tab$ti[k,s] > 0.0 ) {
                if ( ( k == 1 ) || ( neg.ind.rows(tab, k, s) == 0 ) ) {
                    if ( tab$ti[k,s] > vsp ) {
                        sp <- s
                        vsp <- tab$ti[k,s]
                    }
                }
            }
        }    
    }
    return( sp )
}
