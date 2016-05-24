print.llgpcptab <- function( x, ... )
{
###
### This function prints out the augmented tableau at the k-th priority level
###
### Parameters
### x = an object of class 'llgpcptab' that is the modified simplex tableau
### ... = other arguments as they may apply to the generic S3 print function
###
    tab <- x
    fmt.e <- paste( "%", 14, ".", 6, "e", sep="" )
    fmt.s <- paste( "%", 14, "s", sep="" )
    cat( "\n" )
    cat( "Iteration Number: ", tab$iter, "\n" )
    cat( "Priority Level: ", tab$level, "\n" )
###
### create, format and print a data frame for the top stub
###
    cat( "\n" )
    cat( "Top Stub\n" )
    tw.frame <- data.frame( matrix( nrow=tab$levels, ncol=tab$nonbasics ) )
    for ( i in 1:tab$levels ) {
        for ( j in 1:tab$nonbasics ) {
            tw.frame[i,j] <- sprintf( fmt.e, tab$tw[tab$levels-i+1,j] )
        }
    }    
    names( tw.frame ) <- sprintf( fmt.s, tab$col.headings )
    row.names( tw.frame ) <- paste( "P", seq( tab$levels, 1, -1 ), sep="" )
    print( tw.frame )
###
### create, format and print a data frame for the center stub
###
    cat( "\n" )
    cat( "Center Stub\n" )
    te.frame <- data.frame( matrix( nrow=tab$objectives, ncol=tab$nonbasics ) )
    for ( i in 1:tab$objectives ) {
        for ( j in 1:tab$nonbasics ) {
            te.frame[i,j] <- sprintf( fmt.e, tab$te[i,j] )
        }
    }    
    names( te.frame ) <- sprintf( fmt.s, tab$col.headings )
    row.names( te.frame ) <- tab$row.headings
    print( te.frame )
###
### create, format and print a data frame for the left stub
###
    cat( "\n" )
    cat( "Left Stub\n" )
    tu.frame <- data.frame( matrix( nrow=tab$objectives, ncol=tab$levels ) )
    for ( i in 1:tab$objectives ) {
        for ( j in 1:tab$levels ) {
            tu.frame[i,j] <- sprintf( fmt.e, tab$tu[i,tab$levels-j+1] )
        }
    }
    names( tu.frame ) <- sprintf( fmt.s, paste( "P", seq( tab$levels, 1, -1 ), sep="" ) )
    row.names( tu.frame ) <- tab$row.headings
    print( tu.frame )
###
### create, format and print a data frame for the index rows
###
    cat( "\n" )
    cat( "Index Rows\n" )
    ti.frame <- data.frame( matrix( nrow=tab$levels, ncol=tab$nonbasics ) )
    for ( i in 1:tab$levels ) {
        for ( j in 1:tab$nonbasics ) {
            ti.frame[i,j] <- sprintf( fmt.e, tab$ti[i,j] )
        }
    }
    names( ti.frame ) <- sprintf( fmt.s, tab$col.headings )
    row.names( ti.frame ) <- paste( "P", 1:tab$levels, sep="" )
    print( ti.frame )
###
### create, format and print a data frame for the current solution
###
    cat( "\n" )
    cat( "Current Solution\n" )
    tb.frame <- data.frame( matrix( nrow=tab$objectives, ncol=1 ) )
    for ( i in 1:tab$objectives ) {
        tb.frame[i,1] <- sprintf( fmt.e, tab$tb[i] )
    }
    names( tb.frame ) <- sprintf( fmt.s, c( "Value" ) )
    row.names( tb.frame ) <- tab$row.headings
    print( tb.frame )
###
### create, format and print a data frame for the achievement function
###
    cat( "\n" )
    cat( "Achievement Function\n" )
    ta.frame <- data.frame( matrix( nrow=tab$levels, ncol=1 ) )
    for ( i in 1:tab$levels ) {
        ta.frame[i,1] <- sprintf( fmt.e, tab$ta[i] )
    }
    names( ta.frame ) <- sprintf( fmt.s, c( "Value" ) )
    row.names( ta.frame ) <- paste( "P", 1:tab$levels, sep="" )
    print( ta.frame )
###
### print the variable classes
###
    cat( "\nVariable Classes\n" )
    print( tab$variable.classes )
}
