
slope.com <- function( y, x, groups, method="SMA", alpha=0.05, data=NULL, intercept=TRUE, robust=FALSE, V=array( 0, c( 2,2,length(unique(groups)) ) ), 
	group.names=sort(unique(groups)), ci=TRUE, bs=TRUE, slope.test=NULL )
{
    if ( nargs() < 3 )
    {
        stop('Sorry, no can do without three arguments -- Y, X, GROUPS')
    }

#     if ( is.null(data)==FALSE )
#     {
#         attach(data)
#     }

    
    if(!is.null(data))
      stop("'data' argument no longer supported.")
    
    
    dat    <- cbind(y, x)
    g      <- length(group.names)

    # Find sample size, variances for each group:
    n      <- matrix( 0, g, 1 )
    r.factor      <- matrix( 0, g, 1 )
    res.df <- matrix( 0, g, 1 )
    z      <- matrix( 0, g, 3 )
    do.bs  <- bs
    bs     <- matrix( NA, 3, g, dimnames=list(c("slope","lower.CI.lim","upper.CI.lim"),group.names) )
    for (i in 1:g)
    {
        iref   <- ( groups==group.names[i] )
        iref   <- iref & ( is.na(x+y) == FALSE )
        n[i]   <- sum(iref)
        if ( robust )
        {
	   if ( intercept==FALSE )
           {
		stop("Sorry, robust estimation without an intercept term not yet implemented.")
	   }
           else 
           {
                if (n[i]>1)
                    { q     <- pchisq(3,2)
                      S     <- huber.M(dat[iref,])
		      xi     <- S$cov-V[, , i]  
        	      means <- S$loc
                     
		      # get robust.factor for group i (multiplier on variance matrix):
	              r.factor[i] <- robust.factor(dat[iref,],q)[1] 
		    }
                    else if (n[i]==1)
                    { xi <- matrix(0,2,2)
                    r.factor[i] <- 0  } #leave as zero for n[i]=1
           }
           z[i,]     <- c( xi[1,1], xi[2,2], xi[1,2] )

           if (do.bs==TRUE & n[i]>1)
           {
                slopei    <- slope.test(y[iref], x[iref], method=method, alpha=alpha, V=V[,,i], intercept=intercept, robust=TRUE)
                bs[,i]    <- c(slopei$b, slopei$ci)
           }
        }
    	else
    	{
           r.factor[i]<-1
       
           if ( intercept==FALSE )
           {
	  	xi <- t(dat[iref, ]) %*% dat[iref, ] / n[i] - V[, , i]
	   }
           else 
           {
                if (n[i]>1)
                    { xi <- cov(dat[iref, ]) - V[, , i] }
                    else if (n[i]==1)
                    { xi <- matrix(0,2,2) } #leave as zero for n[i]=1
           }

           z[i,]     <- c( xi[1,1], xi[2,2], xi[1,2] )

           if (do.bs==TRUE & n[i]>1)
           {
                slopei    <- slope.test(y[iref], x[iref], method=method, alpha=alpha, V=V[,,i], intercept=intercept,robust=FALSE )
                bs[,i]    <- c(slopei$b, slopei$ci)
           }
        }
    }

    if (intercept==FALSE)
        { res.df <- n-1 }
    else
        { res.df <- n-2 }
#    res.df = res.df + as.numeric( is.null(slope.test)==F ) #to add one to df if slope is given a priori.

#     if ( is.null(data)==FALSE )
#     {
#         detach(data)
#     }

    # Find common slope:
    lambda <- 1 #only actually used for the major axis.
    if (is.null(slope.test))
    {
        res    <- b.com.est( z, n, method, lambda, res.df=res.df )
    }
    else #input slope.test as common slope value.
    {
        if ( (method==1) | (method=='SMA') ) { lambda <- slope.test^2 }
        res <- list(b=slope.test, bi=slope.test, l1=NA, l2=NA, lambda=lambda )
    }
    # Calculate LR:
    dets <- z[,1]*z[,2] - z[,3]^2 #This is l1*l2 under Halt.
    arguments <- list( l1=dets, l2=1, z=z, n=n, method=method, crit=0, lambda=lambda, res.df=res.df, r.factor=r.factor)
    LR     <- lr.b.com(res$b, arguments) 
    # if lambda is being estimated, check endpoint LR values:
    if ( (method==3) | ( method=='lamest' ) )
    {
        res0      <- b.com.est( z, n, 2, lambda=10^-9, res.df ) # to find est when lambda=0
        arguments <- list( l1=dets, l2=1, z=z, n=n, method=method, crit=0, lambda=10^-9, res.df=res.df, r.factor=r.factor)
        LR0       <- lr.b.com(res0$b, arguments) 
        resInf    <- b.com.est( z, n, 2, 10^9, res.df ) # to find est when lambda=inf
        arguments <- list( l1=dets, l2=1, z=z, n=n, method=method, crit=0, lambda=10^9, res.df=res.df, r.factor=r.factor)
        LRinf     <- lr.b.com(resInf$b, arguments) 
        LR        <- min(LR,LR0,LRinf)
        if ( LR==LR0 )    { res <- res0 }
        if ( LR==LRinf )  { res <- resInf }
    }
    
    # Record values for arguments separately
    b      <- res$b
    bi     <- res$bi
    l1     <- res$l1
    l2     <- res$l2
    lambda <- res$lambda

    # Calculate P-value:
    if (is.null(slope.test))
    {
        df <- g - 1 - sum(n<=1)  #don't count any singleton or empty groups in df
    }
    else
    {
        df <- g - sum(n<=1)  #if common slope is given, don't subtract a df for its estimation
    }
    Pvalue <- 1 - pchisq( LR, df )

    # Calculate a CI for common slope, if estimated
    if (is.null(slope.test))
    {
        if ( (method==1) | (method=='SMA') )
        {
           # Getting variance of common slope
           varBs <- ( l2/l1 + l1/l2 + 2)^(-1) *  ( 4*b^2 ) * r.factor
           l12   <- ( l1 + l2 )^2 / ( l1*l2 )
      
        } 
        else
        if ( (method==2) | (method=='MA') )
        {
           varBs <- ( l2/l1 + l1/l2 - 2)^(-1) * ( lambda + b^2 )^2 * r.factor
           l12 <- ( l1 - l2 )^2 / ( l1*l2 )
        }
        if ( (method==3) | (method=='lamest') )
        #Still work to be done to calculate CI for lamest.
        {
           varBs <- NA
        }

        varB <- sum( l12^2 * varBs * res.df )/(sum( l12 * res.df ))^2
    }
    else
    {
        varB <- NA
        ci   <- FALSE
    }

    crit <- qchisq( 1 - alpha, 1 )
    if ( (method==3) | (method=='lamest') )
    {
       ci <- FALSE
    }
    bCI=NA
    if ( ci == TRUE )
    {
       bCI  <- com.ci( b, varB, crit, z, n, l1, l2, method, lambda, res.df, r.factor)
    }
    if (lambda==10^-9) { lambda <- 0 }
    if (lambda==10^9) { lambda  <- Inf }
    list( LR=LR, p=Pvalue, b=b, ci=bCI, varb=varB, lambda=lambda, bs=bs, df=df )
}

