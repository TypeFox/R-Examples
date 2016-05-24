
summary.cvsvd <- function( object, ... ) {
    msep    <- object$msep
    maxrank <- object$maxrank
    
    K         <- nrow( msep )
    rank      <- seq( from=0, to=maxrank, by=1 )
    msep.mean <- apply( msep, 2, mean )
    msep.se   <- apply( msep, 2, sd ) / sqrt( K )
    rank.best <- which.min( msep.mean ) - 1
    rank.1se  <- min( which( msep.mean 
                               <= msep.mean[ rank.best+1 ] 
                                  + msep.se[ rank.best+1 ] ) ) - 1
    
    names( rank.best ) <- NULL
                                  
    list( nfolds=K, maxrank=maxrank, 
          msep.mean=msep.mean, msep.se=msep.se, 
          rank.best=rank.best, rank.1se=rank.1se ) 
}
