
# Bootstrapped CI of SPLS coefficients in active sets

"ci.spls" <-
function( object, coverage=0.95, B=1000,
        plot.it=FALSE, plot.fix='y', plot.var=NA,
        K=object$K, fit=object$fit )
{
    # initialization

    betahat <- object$betahat
    y <- object$y
    A <- object$A
    x <- object$x
    xA <- object$x[,A,drop=FALSE]
    n <- nrow(y)
    p <- ncol(x)
    q <- ncol(y)

    # bootstrap

    betamat <- array( 0, c(length(A),q,B) )

    for ( i in 1:B )
    {
        if ( i%%ceiling(B/10)==0 )
        {
            perc <- round( 10*i/ceiling(B/10) )
            cat( paste(perc,'% completed...\n') )
        }
        nbt <- sample( c(1:n), replace=TRUE )
        ybt <- y[nbt,,drop=FALSE]
        xAbt <- xA[nbt,,drop=FALSE]
        plsfit <- plsr( ybt~xAbt, ncomp=K, method=fit )
        betamat[,,i] <- coef(plsfit)
    }

    # calculate CI

    tailp <- ( 1 - coverage ) / 2
    qt <- function(x) { quantile( x, c(tailp,1-tailp) ) }
    cibeta <- list()
    lbmat <- ubmat <- matrix( 0, p, q )
    for ( i in 1:q )
    {
        cii <- t( apply( betamat[,i,], 1, qt ) )
        cibeta[[i]] <- cii
        rownames(cibeta[[i]]) <- colnames(xA)
        lbmat[A,i] <- cii[,1]
        ubmat[A,i] <- cii[,2]
    }
    names(cibeta) <- colnames(y)

    # CI plot

    if ( plot.it==TRUE )
    {
        if ( plot.fix=='y' )
        {
            if ( is.na(plot.var[1]) ) { plot.var <- c(1:q) }
            k <- 1
            for ( i in plot.var )
            {
                if ( k>1 ) { dev.new() }
                ylimit <- c( min(lbmat[,i]), max(ubmat[,i]) )
                plot( A, betahat[A,i], type='p',
                    xlim=c(1,p), ylim=ylimit,
                    xlab='Predictors',
                    ylab='Coefficient Estimates',
                    main=paste( format(100*coverage,nsmall=0),
                    '% Bootstrapped CI of Coefficients',sep='') )
                for ( j in 1:length(A) )
                { lines( c(A[j],A[j]), c(lbmat[A[j],i],ubmat[A[j],i]) ) }
                abline( h=0, lty=2, col='red' )
                k <- k + 1
            }
        }
        if ( plot.fix=='x' )
        {
            if ( is.na(plot.var[1]) ) { plot.var <- A }
            k <- 1
            for ( i in plot.var )
            {
                if ( k>1 ) { dev.new() }
                ylimit <- c( min(lbmat[i,]), max(ubmat[i,]) )
                plot( c(1:q), betahat[i,], type='p',
                    xlim=c(1,q), ylim=ylimit,
                    xlab='Responses',
                    ylab='Coefficient Estimates',
                    main=paste( format(100*coverage,nsmall=0),
                    '% Bootstrapped CI of Coefficients',sep='') )
                for ( j in 1:q )
                { lines( c(j,j), c(lbmat[i,j],ubmat[i,j]) ) }
                abline( h=0, lty=2, col='red' )
                k <- k + 1
            }
        }
    }

    ci <- list( cibeta=cibeta, betahat=betahat,
                lbmat=lbmat, ubmat=ubmat )
    invisible(ci)
}
