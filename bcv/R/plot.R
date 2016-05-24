
plot.cvsvd <- function( x, errorbars = TRUE, add = FALSE,
                        xlab = "Rank", ylab = "Mean Sq. Prediction Error", 
                        col = "blue", col.errorbars = "gray50", 
                        ... ) {
    msep   <- x$msep
    maxrank <- x$maxrank
    
    K          <- nrow( msep )
    rank       <- seq( from=0, to=maxrank, by=1 )
    msep.mean  <- apply( msep, 2, mean )
    msep.se    <- apply( msep, 2, sd ) / sqrt( K )
    
    if( !add ) {
        if( errorbars ) {
            plot( c(rank-0.2,rank+0.2), msep.mean+c(-msep.se, msep.se), 
                  t='n', xlab=xlab, ylab=ylab, ... )
        } else {
            plot( rank, msep.mean, t='n', xlab=xlab, ylab=ylab, ... )
        }
    }

    lines( rank, msep.mean, col=col, ... )
    
    if( errorbars ) {
        segments( rank-0.2, msep.mean-msep.se, 
                  rank+0.2, msep.mean-msep.se,
                  col=col.errorbars )

        segments( rank, msep.mean-msep.se, 
                  rank, msep.mean+msep.se,
                  col=col.errorbars )

        segments( rank-0.2, msep.mean+msep.se, 
                  rank+0.2, msep.mean+msep.se,
                  col=col.errorbars )
    }

    points( rank, msep.mean, col=col, pch=16, cex=0.6 )
    
    invisible( list( k=rank, msep=msep.mean ) )
}
