
shift.com <- function( y, x, groups, data=NULL, method="SMA", intercept=TRUE, robust=FALSE ,  V=array( 0, c( 2,2,length(unique(groups)) ) ), group.names=sort(unique(groups)))
{
#     if ( is.null(data)==FALSE )
#     {
#         attach(data)
#     }

  
  if(!is.null(data))
    stop("'data' argument no longer supported.")
  
  
    y <- as.matrix(y)
    x <- as.matrix(x)
    dat    <- cbind(y, x)
    groups <- as.matrix(groups)

    nobs <- length(groups)
    g    <- length(group.names)
    inter<- intercept

    res  <- slope.com( y, x, groups, method, intercept=inter, V=V, ci=FALSE, bs=FALSE, robust=robust )
    lr   <- res$lr
    p    <- res$p
    b    <- res$b
    varb <- res$varb

    n        <- matrix( 0, g, 1 )
    varAxis  <- n
    as       <- n
    means    <- matrix( 0, g, 2 )

    if ( (method=="SMA") | method==1 )
    {
        axis       <- y + b*x
        coefV1     <- 1 #The coef of V[1,1,:] in var(axis).
        coefV2     <- b^2 #The coef of V[2,2,:] in var(axis).
        mean.ref   <- 2 #Ref for the column of means to use as coef of var(b)
    }
    if ( (method=="MA") | method==2 )
    {
        axis       <- b*y + x
        coefV1     <- b^2 #The coef of V[1,1,:] in var(axis).
        coefV2     <- 1
        mean.ref   <- 1 #Ref for the column of means to use as coef of var(b)
    }
 
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
                
                if ( (method=="SMA") | method==1 )
                {
			as[i]      <- means[i,1] + b*means[i,2]				
			varAxis[i] <- (S$cov[1,1] + 2*b*S$cov[1,2] + b^2*S$cov[2,2]) * r.factor2
                } 

 		if ( (method=="MA") | method==2 )
    		{
			as[i]      <- b * means[i,1] + means[i,2]				
			varAxis[i] <- (b^2 * S$cov[1,1] + 2*b*S$cov[1,2] + S$cov[2,2]) * r.factor2
 		}
        }
        else
        {
		means[i,1] <- mean( y[iref] ) 
       	      	means[i,2] <- mean( x[iref] )
       	       	as[i]      <- mean( axis[iref] )
       	       	varAxis[i] <- var( axis[iref] )
        }

    }
    varAxis    <- varAxis - coefV1*V[1,1,] - coefV2*V[2,2,]
    varAxis    <- varAxis * (n-1) / (n-2)
    mean.for.b <- means[,mean.ref]

    varAs <- diag( array(varAxis/n) ) + varb*mean.for.b%*%t(mean.for.b)

    varAs[n==1,] <- 0 #For singleton groups
    varAs[,n==1] <- 0
    df     <- g - 1 - sum(n==1)
    L      <- matrix(0,df,g)
    L[,n>1] <- cbind( matrix( 1, df, 1), diag( array( -1, df), nrow=df ) )
    stat  <- t(L%*%as)%*%solve(L%*%varAs%*%t(L), tol=1.0e-050 )%*%(L%*%as)

    pvalue <- 1 - pchisq( stat, df )

#     if ( is.null(data)==FALSE )
#     {
#         detach(data)
#     }

    list( stat=stat, p=pvalue, f.mean=as.vector(as), df=df )

}
