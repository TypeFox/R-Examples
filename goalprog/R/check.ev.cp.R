check.ev.cp <- function( tab, s )
{
###
### This function checks the candidate non-basic variable to determine
### if it violates the complementary pivoting requirement
###
### Parameters
### tab = augmented tableau for LLGP with complementary pivoting
### s = the index of the non-basic variable to enter the solution basis
###
###
### get the corresponding non-basic variable
###
    nonbasic.variable <- tab$col.headings[s]
###
### skip over N and P variables
###
    firstLetter <- substr( nonbasic.variable, 1, 1 )
    if ( firstLetter != "X" )
        return( TRUE )
    ev.class <- get.variable.class( tab, nonbasic.variable )
    for ( i in 1:tab$objectives ) {
        basic.variable <- tab$row.headings[i]
        dv.class <- get.variable.class( tab, basic.variable )
        if ( ev.class == dv.class ) {
            return( FALSE )
        }    
    }        
    return( TRUE )
}
