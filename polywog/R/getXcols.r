##
## Match supplied character string to variable names in data frame, and warn
## about no match or multiple partial matches.  Used in 'predVals', not meant to
## be called directly by users.
##
## ARGUMENTS:
##   vars: character string of supplied names
##   names: character string of names in data
##
## RETURN:
##   Indices of elements in 'names' matching elements in 'vars
##
getXcols <- function(vars, names)
{
    ## Figure out the index of 'names' corresponding to each value in 'vars' (NA
    ## for no match, 0 for multiple partial matches)
    vcols <- charmatch(vars, names)
    if (any(is.na(vcols))) {
        stop("The following variable names were not found in the fitted model: ",
             paste(vars[is.na(vcols)], collapse = ", "))
    } else if (any(vcols == 0)) {
        stop("The following names partially match mutliple variables in the fitted model: ",
             paste(vars[vcols==0], collapse = ", "))
    }

    return(vcols)
}
