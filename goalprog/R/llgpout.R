llgpout <- function( tab, coefficients, targets )
{
###
### This function returns a list of named components with the decision variables,
### positive deviation variables and negative deviation variables.
###
### Parameters
### tab = a list of named components that specifies the modified simplex tableau
### coefficients = a matrix of the coefficients for the linear objective functions
### targets = a vector of the target values for the objective functions
###
### create vectors with the labels for the variables
###
    x.headings <- paste( "X", 1:tab$variables, sep="" )
    neg.headings <- paste( "N", 1:tab$objectives, sep="" )
    pos.headings <- paste( "P", 1:tab$objectives,sep="" )
###
### allocate arrays for the output variables
###
    x <- rep( 0, tab$variables )
    n <- rep( 0, tab$objectives )
    p <- rep( 0, tab$objectives )
###
### loop over the basic variables in the tableau
###
    for ( i in 1:tab$objectives ) {
###
###     loop over the decision variable labels
###
        for ( j in 1:tab$variables ) {
            if ( tab$row.headings[i] == x.headings[j] ) {
                x[j] <- tab$tb[i]
            }
        }
###
###     loop over the negative deviation variable labels
###
        for ( j in 1:tab$objectives ) {
            if ( tab$row.headings[i] == neg.headings[j] ) {
                n[j] <- tab$tb[i]
            }
        }
###
###     loop over the positive deviation variable labels
###
        for ( j in 1: tab$objectives ) {
            if ( tab$row.headings[i] == pos.headings[j] ) {
                p[j] <- tab$tb[i]
            }
        }
    }
    f <- coefficients %*% x
    output <- list( x=x, n=n, p=p, f=f, a=tab$ta, b=targets )
###
### define the llgpout class
###
    class( output ) <- "llgpout"
    return( output )
}
