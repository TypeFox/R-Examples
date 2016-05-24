
"wpls" <-
function( x, y, V, K=ncol(x), type="pls1",
        center.x=TRUE, scale.x=FALSE )
{
    n <- nrow(x)
    p <- ncol(x)
    q <- ncol(y)
    
    # initialization
    
    x1 <- x
    y1 <- y
    
    # standardize x
    
    #if ( center.x )
    #{
    #    x1 <- t(x1) - apply( x1, 2, weighted.mean, diag(V) )
    #    x1 <- t(x1)
    #}
    #if ( scale.x )
    #{
    #    x1 <- t(x1) / sqrt( apply( x1, 2, var ) )
    #    x1 <- t(x1)
    #}    
    
    # main iteration
    
    W <- matrix( 0, p, K )
    T <- matrix( 0, n, K )
    Q <- matrix( 0, q, K )
    P <- matrix( 0, p, K )
    
    for ( k in 1:K )
    {
        # direction vector
        
        #w <- t(x1) %*% V %*% y1
        w <- t(x1) %*% as.matrix( V * y1 )
        #w <- w / sqrt( drop( t(w) %*% w ) )
        w <- w / sqrt( sum(w^2) )
        W[,k] <- w
        
        # latent component
        
        t <- x1 %*% w 
        T[,k] <- t
        
        # standardize t?
        # t <- t - weighted.mean( t, diag(V) )
        # t <- t / sd(t)
        
        # coefficient 
        
        #coef.q <- t(t) %*% V %*% y1 / drop( t(t) %*% V %*% t )
        coef.q <- sum(t*V*y1) / sum( t*V*t )
        Q[,k] <- coef.q
        #coef.p <- t(t) %*% V %*% x1 / drop( t(t) %*% V %*% t )
        coef.p <- t(as.matrix(t*V)) %*% x1 / sum( t*V*t )
        P[,k] <- coef.p
        
        # update
        
        if ( type=='pls1' )
        {
            y1 <- y1 - t %*% coef.q
            x1 <- x1 - t %*% coef.p
        }
        if ( type=='simpls' )
        {
            pj <- w
            pw <- pj %*% solve(t(pj) %*% pj) %*% t(pj)
            x1 <- x1 - x1 %*% pw
            #pw <- w %*% t(w) / drop( t(w) %*% w )
            #x1 <- x1 - x1 %*% pw
        }
    }
    
    list( W=W, T=T, Q=Q, P=P )
}
