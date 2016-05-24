print.llgpout <- function( x, ... )
{
###
### This function does a pretty print of the give optimal solution
###
### Parameters
### x = an object of type 'llgpout' which is the optimal solution
### ... = other arguments as they may apply to the generic S3 print function
###
### create format string for the numerical values
###
    output <- x
    fmt.e <- paste( "%", 14, ".", 6, "e", sep="" )
###
### create format string for a character string
###
    fmt.s <- paste( "%", 14, "s", sep="" )
###
### transform the results in the llgpout object to string vectors
###
    x <- sprintf( fmt.e, output$x )
    n <- sprintf( fmt.e, output$n )
    p <- sprintf( fmt.e, output$p )
    f <- sprintf( fmt.e, output$f )
    b <- sprintf( fmt.e, output$b )
    a <- sprintf( fmt.e, output$a )
###
### get the lengths of these vectors
###
    variables <- length( x )
    objectives <- length( f )
    levels <- length( a )
###
### create the data frames with character strings
###
    solution.frame <- data.frame( matrix( x, nrow=variables, ncol=1 ) )
    names( solution.frame ) <- sprintf( fmt.s, c( "X" ) )
    row.names( solution.frame ) <- paste( "X", 1:variables, sep="" )
    objectives.frame <- data.frame( cbind( f, p, n, b ) )
    names( objectives.frame) <- sprintf( fmt.s,
        c( "Objective", "Over", "Under", "Target" ) )
    row.names( objectives.frame ) <- paste( "G", 1:objectives, sep="" )
    achievement.frame <- data.frame( matrix( a, nrow=levels, ncol=1 ) )
    names( achievement.frame ) <- sprintf( fmt.s, c("A" ) )
    row.names( achievement.frame ) <- paste( "P", 1:levels, sep="" )
###
### print the data frames
###
    cat( "\nDecision variables\n" )
    print( solution.frame )
    cat( "\nSummary of objectives\n" )
    print( objectives.frame )
    cat( "\nAchievement function\n" )
    print( achievement.frame )
}