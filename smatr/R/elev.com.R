
elev.com <- function( y, x, groups, data=NULL, method="SMA", alpha=0.05, robust=FALSE, V=array( 0, c( 2,2,length(unique(groups)) ) ), 
    group.names=sort(unique(groups)) )
{
#     if ( is.null(data)==FALSE )
#     {
#         attach(data)
#     }

    if(!is.null(data)){
      stop("'data' argument no longer supported. Use with() instead.")
    }
  
    x      <- as.matrix( x )
    y      <- as.matrix( y )
    dat    <- cbind(y, x)
    groups <- as.matrix( groups )

    nobs <- length( groups )
    g    <- length( group.names )
    res  <- slope.com( y, x, groups, method=method, V=V, bs=FALSE, ci=FALSE, robust=robust )   
    lr   <- res$lr
    p    <- res$p
    b    <- res$b
    varb <- res$varb

    n      <- matrix( 0, g, 1 )
    varres <- matrix( 0, g, 1 )
    means  <- matrix( 0, g, 2 )
    res    <- y - b*x

    for ( i in 1:g )
    {
        iref       <- ( groups==group.names[i] )
        iref       <- iref & ( is.na(x+y) == FALSE )
        n[i]       <- sum( iref )

        if (robust)       
        {
                q          <- pchisq(3,2)
                S          <- huber.M(dat[iref,])
                r.mean     <- S$loc

                # robust factor for means:
                r.factor2 <- robust.factor(dat[iref,],q)[2]

        	means[i,1] <- r.mean[1] 
        	means[i,2] <- r.mean[2]

                varres[i] <- ((S$cov[1,1] - 2*b*S$cov[1,2] + b^2*S$cov[2,2]) * r.factor2 - V[1,1,i] - b^2*V[2,2,i] )
        }
        else
        {
		means[i,1] <- mean( y[iref] ) 
        	means[i,2] <- mean( x[iref] )
        	varres[i]  <- ( var( res[iref] ) - V[1,1,i] - b^2*V[2,2,i] )
        }
    }

    varres <- varres *( n - 1 )/( n - 2 )

    as     <- means[,1] - b*means[,2]
    names(as) <- group.names
    varas  <- diag( array( varres/n ) ) + varb * means[,2]%*%t( means[,2] )

    varas[n==1,] <- 0 #For singleton groups
    varas[,n==1] <- 0
    df     <- g - 1 - sum(n==1)
    L      <- matrix(0,df,g)
    L[,n>1] <- cbind( matrix( 1, df, 1), diag( array( -1, df), nrow=df ) )
    
    stat   <- t(L%*%as)%*%( solve( L%*%varas%*%t(L) ) )%*%(L%*%as)
    pvalue <- 1 - pchisq( stat, df ) 
    sinv=matrix(0,g,g)
    sinv[n>1,n>1]   <- solve( varas[n>1,n>1] )
    a      <- (matrix(1,1,g)%*%sinv%*%as)/ sum( sum( sinv) )
    vara   <- 1 / sum( sum( sinv ) )
    crit   <- qchisq( 1 - alpha, 1 )
    crits  <- qt( 1 - alpha/2, n-2 )
    a.ci   <- c( a - sqrt( crit*vara), a + sqrt( crit*vara ) )
    as.ci  <- as + cbind(0, - crits*sqrt(diag(varas)), crits*sqrt(diag(varas)) )
	
    dimnames(as.ci)[[1]] = group.names
    dimnames(as.ci)[[2]] <- c("elevation","lower CI limit","upper CI limit")
# 
#     if ( is.null(data)==FALSE )
#     {
#         detach(data)
#     }

    list( stat=stat, p=pvalue , a=a, ci=a.ci, as = as.ci, df=df )
}
