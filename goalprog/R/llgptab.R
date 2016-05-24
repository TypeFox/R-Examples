llgptab <- function( coefficients, targets, achievements )
{
###
### This function returns a list of named components that specify
### the initial modified simplex tableau for the lexicographical
### linear programming model
###
### Parameters
### coefficients = a matrix with the coefficients of the linear objective functions
### targets = a vector of target values for the objective functions
### achievements = a data frame with the weights of the deviation variables for each
### objective function along with the corresponding priority level
###
###        
### get the highest priority level from the achievements data frame
###
    levels <- max( achievements$priority )
###
### ensure that there is at least one objective that is an absolute objective
###
#    if ( min( achievements$priority ) != 1 )
#        stop( "there is no absolute objective in the achievements" )
###
### get the number of goals from the achievements data frame
###
    goals <- nrow( achievements )
###
### get the number of variables  and objectivesfrom the coefficients matrix
###
    variables <- ncol( coefficients )
    objectives <- nrow( coefficients )
###
### compute the number of non basic variables
###
    nonbasics <- variables + objectives
###
### ensure dimensional consistency among the arguments
###
    if ( nrow( coefficients ) != objectives )
        stop( "coefficients and achievements do not have the same number of objectives " )
    if ( length( targets ) != objectives )
        stop( "achievements and targets do not have the same number of objectives" )
###
### create matrix of coefficients for the tableau
###
    te <- cbind( coefficients, diag( -1.0, objectives ) )
###
### create the matrices for the top and left stubs
###
    tw <- matrix( 0, nrow=levels, ncol=nonbasics )
    tu <- matrix( 0, nrow=objectives, ncol=levels )
###
### transfer achievements to the stubs
###
    for ( goal in 1:goals ) {
        i <- achievements[goal,"objective"]
        k <- achievements[goal,"priority"]
        w <- achievements[goal,"p"]
        u <- achievements[goal,"n"]
        s <- variables + i
        tw[k,s] <- w
        tu[i,k] <- u
    }
###
### create row headings
###
    col.headings <- c( paste( "X", 1:variables, sep="" ), paste( "P", 1:objectives,sep="" ) )
    row.headings <- paste( "N", 1:objectives, sep="" )
###
### initialize the target and achievement vectors of the tableau
###
    tb <- targets
    ta <- rep( 0, levels )
###
### create the matrix of the index rows
###
    ti <- matrix( 0, nrow=levels, ncol=nonbasics )
###
### create the list with the tableau
###
    tab <- list( iter=0, variables=variables, levels=levels, objectives=objectives,
                 nonbasics=nonbasics, level=0, te=te, tb=tb, tw=tw, tu=tu, ti=ti, ta=ta,
                 row.headings=row.headings, col.headings=col.headings )
###
### define the llgptab class
###
    class( tab ) <- "llgptab"
###
### return the tableau
###
    return( tab )
}
