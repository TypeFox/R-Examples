calc.ta.k <- function( tab, k )
{
###
### This function calculates the achievement function for the highest k priority levels
###
### Parameters
### tab = a list of named components that specifies the modified simplex tableau
### k = an integer priority level
###
    tab$ta <- t( tab$tu ) %*% tab$tb
    return( tab )
}
