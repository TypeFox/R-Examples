
# Heatmap-type MSPE plot for s*K

"cv.spls" <-
function( x, y, fold=10, K, eta, kappa=0.5,
        select="pls2", fit="simpls",
        scale.x=TRUE, scale.y=FALSE, plot.it=TRUE )
{
    # initialization
    
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    ip <- c(1:p)
    y <- as.matrix(y)
    q <- ncol(y)
    
    type <- correctp( x, y, eta, K, kappa, select, fit )
    eta <- type$eta
    K <- type$K
    kappa <- type$kappa
    select <- type$select
    fit <- type$fit
    
    # CV MSPE estimation
    
    # data partition
    
    foldi <- split( sample(1:n), rep(1:fold, length = n) )
    
    # heatmap-type plot for multi s & K
        
    mspemat <- matrix( 0, length(eta), length(K) )
    
    for ( i in 1:length(eta) )
    {
        # eta
        
        cat( paste('eta =',eta[i],'\n') )
        
        mspemati <- matrix( 0, fold, length(K) )
        for ( j in 1:fold )
        {
            # fold
            
            #print( paste('fold ',j,sep='') )
            
            omit <- foldi[[j]]
            object <- spls( x[-omit,,drop=FALSE], y[-omit,,drop=FALSE], eta=eta[i], kappa=kappa,
                        K=max(K), select=select, fit=fit,
                        scale.x=scale.x, scale.y=scale.y, trace=FALSE )
            newx <- x[omit,,drop=FALSE]
            newx <- scale( newx, object$meanx, object$normx )
            betamat <- object$betamat
            for (k in K)
            {                    
                # K
                # calculate MSPE
                
                pred <- newx %*% betamat[[k]] + matrix(1,nrow(newx),1) %*% object$mu
                mspemati[j,(k-min(K)+1)] <- mean( apply( (y[omit,]-pred)^2, 2, mean) )
            }
        }
    mspemat[i,] <- apply(mspemati,2,mean)
    }
    
    # find optimal eta & K

    minpmse <- min(mspemat)
    rownames(mspemat) <- eta
    colnames(mspemat) <- K
    mspecol <- apply( mspemat, 2, min)
    msperow <- apply( mspemat, 1, min)
    K.opt <- min( K[mspecol==minpmse] )
    eta.opt <- max( eta[msperow==minpmse] )
    
    cat( paste('\nOptimal parameters: eta = ',eta.opt,', ',sep='') )
    cat( paste('K = ',K.opt,'\n',sep='') )
    
    # plot heatmap & return values
    
    if ( plot.it )
    { heatmap.spls( t(mspemat), xlab='K', ylab='eta', main='CV MSPE Plot', coln=16, as='n' ) }
    rownames(mspemat) <- paste('eta=',eta)
    colnames(mspemat) <- paste('K =',K)
    
    cv <- list( mspemat=mspemat, eta.opt=eta.opt, K.opt=K.opt )
    invisible(cv)
}
