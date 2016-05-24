calc.ti.k <- function( tab, k )
{
###
### This function calculates the index rows for the highest k priority levels
###
### Parameters
### tab = a list of named components that specifies the modified simplex tableau
### k = an integer priority level
###
    tab$ti <- t( tab$tu ) %*% tab$te - tab$tw
    return( tab )
}
