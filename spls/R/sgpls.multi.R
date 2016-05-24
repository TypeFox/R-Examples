
"sgpls.multi" <-
function( x, y, K, eta, scale.x=TRUE,
        eps=1e-5, denom.eps=1e-20, zero.eps=1e-5, maxstep=100,
        br=TRUE, ftype='iden' )
{    
    # initialization
    
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    J <- length(unique(y)) - 1
    ip <- c(1:p)
    p <- ncol(x)
    y.org <- as.matrix(y)+1
    q <- ncol(y.org)
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
    
    ymat <- matrix( 0, n, (J+1) )
    ymat[ cbind( seq(n), y.org ) ] <- 1
    y.base <- ymat[ , 1, drop=FALSE ]
    ymat <- ymat[ , -1 ]
    
    # main iteration
    
    beta1hat <- matrix( 0, p, J )
    beta0hat <- matrix( 0, J, 1 ) 
    beta1hat.old <- beta1hat + 1000
    
    A.bygroup <- list()
    W.bygroup <- list()
    
    re <- 100   # relative error
    min.re <- 1000
    nstep <- 0  # number of iterations
    nstep.min <- 0
    
    groupj <- 1
    
    while ( re > eps & nstep < maxstep )
    {                        
        # estimated success probability & weight
        
        if ( nstep > 0 )
        {
            exp.xb <- exp( t( beta0hat %*% one ) + x0 %*% beta1hat )            
            p0 <- exp.xb[ , groupj ] / ( 1 + apply( exp.xb, 1, sum ) )
            p0[ exp.xb[ , groupj ]==Inf ] <- 1 - zero.eps
            p0[ p0 < zero.eps ] <- zero.eps
            p0[ p0 > (1 - zero.eps) ] <- 1 - zero.eps            
            V <- as.vector( p0 * ( 1 - p0 ) )
            
        } else
        {
            p0 <- 2 * ( ymat[,groupj] + 0.5 ) / ( J + 3 )
            V <- as.vector( p0 * ( 1 - p0 ) )
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
        
        if ( nstep==0 )
        {
            y0 <- matrix( beta0hat[groupj] + x0 %*% beta1hat[,groupj] +
                ( ymat[,groupj] - p0 ) / V )                     
        } else
        {
            Hstar <- (J+1) * H
            Hw <- H - p0 * Hstar
            V <- V * ( Hstar*br/2 + 1 )            
            y0 <- matrix( beta0hat[groupj] + x0 %*% beta1hat[,groupj] +
                ( ymat[,groupj] + H*br/2 - ( Hstar*br/2 + 1 ) * p0 ) / V )
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
        
        A.bygroup[[groupj]] <- A
        W.bygroup[[groupj]] <- W
        
        # update coefficient
        
        beta1hat.old[,groupj] <- beta1hat[,groupj]
        beta1hat[,groupj] <- 0
        beta1hat[A,groupj] <- W %*% solve( t(P) %*% W ) %*% t(Q)      
        beta0hat[groupj] <- weighted.mean( ( y0 - T %*% t(Q) ), sqrt(V) )
        
        # update convergence criterion
        
        re <- mean( abs( beta1hat - beta1hat.old ) ) /
            mean( abs(beta1hat.old) + denom.eps )
        
        if ( re < min.re & nstep > 1 )
        {
            ## record minimum values
            
            min.re <- re
            nstep.min <- nstep
            beta1hat.min <- beta1hat
            beta0hat.min <- beta0hat
            A.bygroup.min <- A.bygroup
            W.bygroup.min <- W.bygroup
        }
        
        # move around group index
        
        if ( groupj < J )
        {
            groupj <- groupj + 1
        } else if ( groupj == J )
        {
            groupj <- 1
            nstep <- nstep + 1
        }
    }
    
    # check convergence
    
    if ( re > eps )
    {
        if ( nstep.min > 0 )
        {
            beta1hat <- beta1hat.min
            beta0hat <- beta0hat.min
            A.bygroup <- A.bygroup.min
            W.bygroup <- W.bygroup.min
        }
    }
    
    # final estimate: active set A
    
    A.all <- sort( unique( unlist( A.bygroup ) ) )
    names(A.bygroup) <- paste( '0.vs.',1:J,sep='' )
    
    # final estimate: coefficients 

    betahat <- rbind( t(beta0hat), beta1hat )
            
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
        betahat=betahat, A.all=A.all, A.bygroup=A.bygroup,
        W.bygroup=W.bygroup, mu=mu, sigma=sigma )
    class(object) <- "sgpls"
    object
}
