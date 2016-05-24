get.variable.class <- function ( tab, variable )
{
###
### This function returns the variable complementarity class for the given variable
###
### Parameters
### tab = a list of named components for the augmented modified simplex tableau
### variable = a character string for the name of a variable
###
    variable.class <- 0
    for ( i in 1:tab$variables ) {
        if ( variable == tab$variable.classes[i,"variable"] ) {
            variable.class <- tab$variable.classes[i,"class"]
        }
    }
    return ( variable.class )
}
