
# fit SPLS direction vector

"spls.dv" <-
function( Z, eta, kappa, eps, maxstep )
{        
    # initialization
    
    p <- nrow(Z)
    q <- ncol(Z)
    Znorm1 <- median( abs(Z) )
    Z <- Z / Znorm1
    
    # main iterations
    
    if ( q==1 )
    {
        # if univariate response, then just soft thresholding
        
        c <- ust( Z, eta )
    }
    
    if ( q>1 )
    {
        # if multivariate response
        
        M <- Z %*% t(Z)
        dis <- 10
        i <- 1
        
        # main iteration: optimize c and a iteratively
        
        # use svd solution if kappa==0.5
        
        if ( kappa==0.5 )
        {        
            # initial value for a & c (outside the unit circle)
            
            c <- matrix( 10, p, 1 )
            c.old <- c
                
            while ( dis>eps & i<=maxstep )
            {
                # optimize a for fixed c
                            
                mcsvd <- svd( M%*%c ) 
                a <- mcsvd$u %*% t(mcsvd$v)
                
                # optimize c for fixed a
                # soft thresholding ( assuming lambda2 -> Inf )
                
                c <- ust( M%*%a, eta )
                
                # calculate discrepancy between a & c
                
                dis <- max( abs( c - c.old ) )                
                c.old <- c
                i <- i + 1
            }
        }
        
        # solve equation if 0<kappa<0.5
        
        if ( kappa>0 & kappa<0.5 )
        {        
            kappa2 <- ( 1 - kappa ) / ( 1 - 2*kappa )
            
            # initial value for c (outside the unit circle)
            
            c <- matrix( 10, p, 1 )
            c.old <- c 
            
            # define function for Lagrange part
            
            h <- function(lambda)
            {
                alpha <- solve( M + lambda*diag(p) ) %*% M %*% c
                obj <- t(alpha) %*% alpha - 1/kappa2^2
                return(obj)
            }
            
            # control size of M & c if too small
            
            if ( h(eps) * h(1e+30) > 0 )
            { while ( h(eps) <= 1e+5 ) { M <- 2*M; c <- 2*c; } }
            
            while ( dis>eps & i<=maxstep )
            {
                # control size of M & c if too small
                
                if ( h(eps) * h(1e+30) > 0 )
                { while( h(eps) <= 1e+5 ) { M <- 2*M; c <- 2*c; } }
                
                # optimize a for fixed c

                lambdas <- uniroot( h, c( eps, 1e+30 ) )$root
                a <- kappa2 * solve( M + lambdas * diag(p) ) %*% M %*% c
                
                # optimize c for fixed a                
                # soft thresholding ( assuming lambda2 -> Inf )
                
                c <- ust( M%*%a, eta )
                
                # calculate discrepancy between a & c
                
                dis <- max( abs( c - c.old ) )                
                c.old <- c
                i <- i + 1         
            }
        }
    }    
    
    return(c)
}
