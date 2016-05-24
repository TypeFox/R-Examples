
print.cvsvd <- function( x, digits=max(3, getOption("digits") - 3), ... ) {
    cat( "\nCall:\n", deparse( x$call ), "\n\n", sep = "" )
    
    msep   <- x$msep
    maxrank <- x$maxrank
    
    K           <- nrow( msep )
    rank        <- seq( from=0, to=maxrank, by=1 )
    msep.mean   <- apply( msep, 2, mean )
    msep.se     <- apply( msep, 2, sd ) / sqrt( K )
    min.rank    <- which.min( msep.mean ) - 1
    min.rank.se <- min( which( msep.mean 
                               <= msep.mean[ min.rank+1 ] 
                                  + msep.se[ min.rank+1 ] ) ) - 1
    
    rank.fmt   <- format( rank )
    rank.width <- max( nchar( rank.fmt[ 1 ] ), nchar( "Rank" ) + 1 )
    
    mean.fmt   <- format( msep.mean, digits=digits )
    mean.width <- max( nchar( mean.fmt[ 1 ] ), nchar( "MSEP" ) + 1 )
    
    se.fmt     <- format( msep.se, digits=digits )
    se.width   <- max( nchar( se.fmt[ 1 ] ), nchar( "SE" ) + 1 )
    
    fmt <- paste( " %", rank.width, "s", 
                  "  %", mean.width, "s",
                  "  %", se.width, "s", sep="" )
    cat( sprintf( fmt, "Rank", "MSEP", "SE" ), "\n", sep="" )
    cat( " " )
    cat( do.call( paste, 
         c( as.list( rep("-", 
                rank.width + 2 + mean.width + 2 + se.width) ), 
            sep='' ) ) )
    cat( "\n" )
    
    for( i in seq_len( maxrank+1 ) ) {
        rank <- i - 1
        mean <- msep.mean[ i ]
        se   <- msep.se[ i ]
        cat( sprintf( fmt, rank.fmt[ i ], mean.fmt[ i ], se.fmt[ i ] ) )
        if ( rank == min.rank && rank == min.rank.se ) {
            cat( " *+ " )
        } else if ( rank == min.rank ) {
            cat( " * " )
        } else if ( rank == min.rank.se ) {
            cat( " + " )
        }
        cat( "\n" )
    }
    
    cat( "\n" )
    
    invisible( x )
}
