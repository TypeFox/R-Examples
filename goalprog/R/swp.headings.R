swp.headings <- function( tab, nr, nc )
{
###
### This function swaps row and column headings
###
### Parameters
### tab = a list of named components that specifies the modified simplex tableau
### nr = an integer row subscript
### nc = an integer column subscript
###
    temp <- tab$col.headings[nc]
    tab$col.headings[nc] <- tab$row.headings[nr]
    tab$row.headings[nr] <- temp
    return( tab )
}
