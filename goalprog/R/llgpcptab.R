llgpcptab <- function( coefficients, targets, achievements, variable.classes )
{
###
### This function returns a list of named components that specify
### the initial modified simplex tableau for the lexicographical
### linear programming model with complementary pivoting
###
### Parameters
### coefficients = a matrix with the coefficients of the linear objective functions
### targets = a vector of target values for the objective functions
### achievements = a data frame with the weights of the deviation variables for each
### objective function along with the corresponding priority level
### variable.classes = a data frame that defines complementarity classes for all 
### the decision variables
###
###
### get the basic llgptab object
###
    tab <- llgptab( coefficients, targets, achievements )
###
### make sure that all variables are defined in the variable.classes frame
###
    variables <- paste( "X", seq(1,tab$variables, 1), sep="" )
    if ( nrow( variable.classes ) != tab$variables )
        stop( "data frame variables.classes non-conformable with the tableau" )
    for ( i in 1:tab$variables ) {
        if ( variables[i] != variable.classes[i,"variable"] )
            stop( paste( "variable", variables[i], " not defined in variable.classes" ) )
    }        
###
### create the list with the augmented tableau
###
    cptab <- list( iter=0, variables=tab$variables, levels=tab$levels, 
        objectives=tab$objectives, nonbasics=tab$nonbasics, level=0, 
        te=tab$te, tb=tab$tb, tw=tab$tw, tu=tab$tu, ti=tab$ti, ta=tab$ta,
        row.headings=tab$row.headings, col.headings=tab$col.headings,
        variable.classes=variable.classes)
###
### define the llgptab class
###
    class( cptab ) <- "llgpcptab"
###
### return the tableau
###
    return( cptab )
}
