frobenius.norm <- function( x )
{
###
### This function is a wrapper function that uses the entrywise.norm function
### with p = 2 to obtain the fronbenius norm
###
### argument
### x = a numeric matrix or vector
###
    return( entrywise.norm( x, 2 ) )
}
