
# coefficient path plot

"plot.spls" <-
function( x, yvar=c(1:ncol(x$y)), ... )
{
    # initialization

    betamat <- x$betamat
    K <- x$K
    eta <- x$eta
    p <- ncol(x$x)

    # coefficient plot

    Ks <- c(1:K)
    for ( i in 1:length(yvar) )
    {
        if ( i>1 ) { dev.new() }
        betamatq <- c()
        for ( j in Ks )
        {
            betamatq <- cbind( betamatq, betamat[[j]][,yvar[i]] )
        }
        ylimit <- c( min(betamatq), max(betamatq) )
        main.name <- paste('Coefficient Path Plot (eta=',eta,')',sep='')
        plot( Ks, Ks, xlim=c(1,K), ylim=ylimit, type='n',
            xlab='K', ylab='Coefficient Estimates', main=main.name, ... )
        if ( length(Ks)>1 )
        {
            for (j in 1:p) { lines( Ks, betamatq[j,], col=j ) }
        } else
        {
            points( rep( 1, length(betamatq) ), betamatq )
        }
        abline( h=0, lty=2, col='red' )
    }
}
