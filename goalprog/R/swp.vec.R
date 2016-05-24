swp.vec <- function (tab, nr, nc)
{
###
### This function swaps row and column vectors between the top and left stub
###
### Parameters
### tab = a list of named components that specifies the modified simplex tableau
### nr = an integer row index
### nr = an integer column index
###
    for ( k in 1:tab$levels ) {
        temp <- tab$tu[nr,k]
        tab$tu[nr,k] <- tab$tw[k,nc]
        tab$tw[k,nc] <- temp
    }
    return( tab )
}
