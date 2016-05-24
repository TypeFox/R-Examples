piv.llgp <- function( tab, nevc, ndvr, verbose )
{
###
### This function performs a simplex pivot to change the basis variables
###
### Parameters:
### tab = a list of named components that specifies the modified simplex tableau
### nevc = the integer column index of the entering variable
### ndvr = the integer row index of the departing variable
### verbose = a logical variable which if true prints the basis change
###
    objectives <- tab$objectives
    nonbasics <- tab$nonbasics
###
### print the basis change
###
    if ( verbose ) {
        cat( "\n" )
        cat( paste( "Iteration:", tab$iter, 
                      ", Entering Variable:", tab$col.headings[nevc],
                      ", Departing Variable:", tab$row.headings[ndvr] ) )
        cat( "\n" )              
    }                  
###
### swap headings
###
    tab <- swp.headings( tab, ndvr, nevc )
###
### swap vectors
###
    tab <- swp.vec( tab, ndvr, nevc )
###
### compute new elements matrix
###
    piv <- tab$te[ndvr,nevc]
    pib <- tab$tb[ndvr]
    for ( i in 1:objectives ) {
        if ( i != ndvr ) {
            pix <- tab$te[i,nevc] / piv
            tab$tb[i] <- ( tab$tb[i] - pix * pib )
            for ( s in 1:nonbasics ) {
                if ( s != nevc ) {
                    tab$te[i,s] <- ( tab$te[i,s] - tab$te[ndvr,s] * pix )
                }
            }
        }
    }
    for ( s in 1:nonbasics ) {
        tab$te[ndvr,s] <- ( tab$te[ndvr,s] / piv )
    }
    for ( i in 1:objectives ) {
        tab$te[i,nevc] <- ( - tab$te[i,nevc] / piv )
    }
    tab$tb[ndvr] <- tab$tb[ndvr] / piv
    tab$te[ndvr,nevc] <- ( 1 / piv )
###
### return the updated tableau
###
    return( tab )
}
