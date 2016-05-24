chkVariation <- function(x)
## Checks the presence of variation in vector 'x'. Returns FALSE if there is no
## variation.
{
    return(length(unique(x)) > 1)
}