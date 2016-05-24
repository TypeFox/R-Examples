check.tb <- function( tab )
{
###
### This function checks for negative target values and repairs the tableau
###
### Parameter
### tab = a list of named components that specifies the modified simplex tableau
###
    for ( i in 1:tab$objectives ) {
        if ( tab$tb[i] < 0.0 ) {
            for ( j in 1:tab$variables ) {
                tab$te[i,j] <- -tab$te[i,j]
            }
            tab <- swp.headings( tab, i, i + tab$variables )
            tab <- swp.vec( tab, i, i + tab$variables )
        }
    }
    return( tab )
}
