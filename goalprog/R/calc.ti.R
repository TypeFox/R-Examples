calc.ti <- function( tab, k )
{
###
### Compute the k-th index row
###
### Parameters
### tab = a list of named components that specify the modified simplex tableau
### k = an integer priority level
###
    for ( s in 1:tab$nonbasics ) {
        sum <- -tab$tw[k,s]
        for ( i in 1:tab$objectives ) {
            sum <- sum + tab$te[i,s] * tab$tu[i,k]
        }
        tab$ti[k,s] <- sum
    }
}
