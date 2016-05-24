
# coefficient plot

"coefplot.spls" <-
function( object, nwin=c(2,2), xvar=c(1:length(object$A)), ylimit=NA )
{
    q <- ncol(object$y)

    if (q==1)
    {
        cat("Sorry! coef.spls is designed only for the multivariate responses.\n")
    } else
    {
        ndiv <- nwin[1]*nwin[2]
        A <- object$A
        betahatA <- object$betahat[A,,drop=FALSE]
        x <- object$x
        xA <- x[,A,drop=FALSE]
        xAname <- colnames(xA)

        if ( is.na(ylimit[1]) )
        {
            ylimit <- c( min(betahatA), max(betahatA) )
        }
        for ( i in xvar )
        {
            if ( i==1 ) { split.screen( nwin ); }
            if ( i>1 & i%%ndiv==1 ) { dev.new(); split.screen( nwin ); }
            if ( i%%ndiv>0 ) { screen( (i%%ndiv) ) } else { screen( ndiv ) }
            plot( c(1:q), betahatA[i,], type='l', ylim=ylimit,
                xlab='Responses', ylab='Coefficient Estimates', main=xAname[i] )
            abline( h=0, lty=2, col='red' )
            if ( i%%ndiv==0 ) { close.screen( all.screens = TRUE ) }
        }
        if ( i%%ndiv>0 ) { close.screen( all.screens = TRUE ) }
    }
}
