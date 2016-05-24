
"sgpls.binary" <-
function( x, y, K, eta, scale.x=TRUE,
        eps=1e-5, denom.eps=1e-20, zero.eps=1e-5, maxstep=100,
        br=TRUE, ftype='iden' )
{        
    # initialization
    
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    ip <- c(1:p)
    y <- as.matrix(y)
    q <- ncol(y)
    one <- matrix(1,1,n)    
        
    # stadardize predictor matrix
    
    mu <- apply( x, 2, mean )    
    x0 <- scale( x, mu, FALSE )    
    
    if ( scale.x )
    {
        sigma <- apply( x, 2, sd )
        x0 <- scale( x0, FALSE, sigma )
    } else
    {
        sigma <- rep( 1, ncol(x) )
        x0 <- x0
    }
        
    # main iteration
    
    beta1hat <- matrix( 0, p, q )
    beta1hat.old <- beta1hat + 1000
    beta0hat <- 0
    re <- 100   # relative error
    min.re <- 1000
    nstep <- 0  # number of iterations
    nstep.min <- 0
    
    while ( re > eps & nstep < maxstep )
    {                                
        # estimated success probability & weight
        
        if ( nstep == 0 )
        {
            p0 <- ( y + 0.5 ) / 2
            V <- as.vector( p0 * (1-p0) )
            A <- c(1:p)
        } else
        {
            exp.xb <- exp( beta0hat + x0 %*% beta1hat )
            p0 <- exp.xb / ( 1 + exp.xb )
            p0[ exp.xb==Inf ] <- 1 - zero.eps
            p0[ p0 < zero.eps ] <- zero.eps
            p0[ p0 > (1 - zero.eps) ] <- 1 - zero.eps
            V <- as.vector( p0 * (1-p0) )
        }
        
        # hat matrix

        switch( ftype,
            hat = {
                H <- hat( sweep( cbind( rep(1,n), x0 ), 1, sqrt(V), "*" ),
                    intercept=FALSE )
            },
            iden = {
                H <- rep( 1, n )
            }
        )        
    
        # working response
        
        if ( nstep==0 )
        {
            y0 <- beta0hat + x0 %*% beta1hat + ( y - p0 ) / V
        } else
        {
            V <- V * ( H*br + 1 )
            y0 <- beta0hat + x0 %*% beta1hat +
                ( y + H*br/2 - (H*br+1) * p0 ) / V
        }
                       
        # PLS step
        
        y1 <- y0
        y1 <- y1 - mean(y1)
        x1 <- x0
        
        A.old <- c()
        
        for (k in 1:K)
        {        
            # define Z
        
            Z <- t(x1) %*% as.matrix( V * y1 )
            
            # fit direction vector
            
            Znorm1 <- median(abs(Z))
            Z <- Z / Znorm1
            what <- ust(Z, eta)
            
            # construct A
            
            A <- sort( unique( c( A.old, ip[ what!=0 ] ) ) )
            
            # fit pls with predictors in A
            
            x0A <- x0[ , A, drop=FALSE ]           
            plsfit <- wpls( x0A, y0, V, K=min(k,length(A)), type="pls1",
                center.x=FALSE, scale.x=FALSE )
            
            # update
            
            A.old <- A            
            y1 <- y0 - plsfit$T %*% t(plsfit$Q)
            x1 <- x0
            x1[,A] <- x0[,A] - plsfit$T %*% t(plsfit$P)
        }
        
        x0A <- x0[ , A, drop=FALSE ]     
        plsfit <- wpls( x0A, y0, V, K=min(K,length(A)), type="pls1",
                center.x=FALSE, scale.x=FALSE )
        
        W <- plsfit$W
        T <- plsfit$T
        P <- plsfit$P
        Q <- plsfit$Q        
        
        # update coefficient
                        
        beta1hat.old <- beta1hat
        beta1hat <- matrix( 0, p, q )
        beta1hat[A,] <- W %*% solve( t(P) %*% W ) %*% t(Q)
        beta0hat <- weighted.mean( ( y0 - T%*%t(Q) ), sqrt(V) )
                
        # update convergence criterion
        
        re <- mean( abs( beta1hat - beta1hat.old ) ) / mean( abs(beta1hat.old) + denom.eps )
        nstep <- nstep + 1        
        
        if ( re < min.re & nstep > 1 )
        {
            ## record minimum values
            
            min.re <- re
            nstep.min <- nstep
            beta1hat.min <- beta1hat
            beta0hat.min <- beta0hat
            A.min <- A
            W.min <- W
        }
    }
    
    # check convergence
    
    if ( re > eps )
    {
        if ( nstep.min > 0 )
        {
            converged <- FALSE
            beta1hat <- beta1hat.min
            beta0hat <- beta0hat.min
            A <- A.min
            W <- W.min
        }
    }      
    
    # final estimate: coefficients 

    betahat <- matrix( c( beta0hat, beta1hat ) )
            
    # return objects
            
    if ( !is.null(colnames(x)) )
    {
        rownames(betahat) <- 1:nrow(betahat)
        rownames(betahat)[1] <- 'intercept'
        rownames(betahat)[2:nrow(betahat)] <- colnames(x)
    } else
    {
        rownames(betahat) <- c( 0, paste("x",1:p,sep="") )
        rownames(betahat)[1] <- 'intercept'
    }
    
    object <- list( x=x, y=y, x0=x0, eta=eta, K=K,
        betahat=betahat, A=A, W=W, mu=mu, sigma=sigma )
    class(object) <- "sgpls"
    object
}
