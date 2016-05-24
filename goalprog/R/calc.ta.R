calc.ta <- function( tab, k )
{
###
### This function calculates the achievement function for the k-th priority level
###
### Parameters
### tab = a list of named components that specifies the modified simplex tableau
### k = an integer priority level
###
    sum <- 0
    for ( i in 1:tab$objectives ) {
        sum <- sum + tab$tb[i] * tab$tu[i,k]
    }
    tab$ta[k] <- sum
}
