
slope.test <- function( y, x, test.value=1, data=NULL, method="SMA", alpha=0.05, V=matrix(0,2,2), intercept=TRUE, robust=FALSE )
{

    if ( nargs() < 2 ) 
    {
        stop('Sorry, no can do without two arguments -- Y, X')
    }

#     if ( is.null(data)==FALSE )
#     {
#         attach(data)
#     }

    
    if(!is.null(data))
      stop("'data' argument no longer supported.")
    
    
    iref <- ( is.na(x+y) == FALSE ) #to remove NA cases
    n    <- sum(iref)

    if ( intercept==FALSE )
    {
        resDF <- n - 1 
    }
    else 
    {
        resDF <- n - 2
    }

    fCrit <- qf( 1-alpha, 1, resDF )

    dat <- cbind(y[iref], x[iref])

    if ( robust )
    {
	    if( intercept )
	    {
		# get robust mean/var matrix:
		q     <- pchisq(3,2)
		S     <- huber.M(dat)
		means <- S$loc
		vr    <- ( S$cov - V) *(n-1)

	        r.factor <- robust.factor(dat,q)[1]
	    }
	    else
	    {
		stop("Sorry, robust estimation without an intercept term not yet implemented.")
	    }
    }
    else
    {
          r.factor <- 1
	  if ( intercept )
     	  {
		vr <- ( cov(dat) - V )*(n-1)
    	  }
    	  else
	  {
		vr <- t(dat)%*%dat - V*n
     	  }
    }

    r <- vr[1,2]/sqrt( vr[1,1]*vr[2,2] )

    bCI     <- matrix( NA, 1, 2 )
    varTest <- matrix( 0, 2, 2 )

    if ( (method==0) | (method=='OLS') )
    {
        b            <- vr[1,2]/vr[2,2]
        varRes       <- ( vr[1,1] - 2*b*vr[1,2] + b^2*vr[2,2] )/resDF
        varB         <- varRes/vr[2,2] * r.factor
        bCI[1,1]     <- b - sqrt(varB)*sqrt(fCrit)
        bCI[1,2]     <- b + sqrt(varB)*sqrt(fCrit)
        varTest[1,1] <- vr[1,1] - 2*test.value*vr[1,2] + test.value^2*vr[2,2]
        varTest[1,2] <- vr[1,2] - test.value*vr[2,2]
        varTest[2,2] <- vr[2,2]
    }
    else if ( (method==1) | (method=='SMA') )
    {
        b            <- sign(vr[1,2])*sqrt(vr[1,1]/vr[2,2])
        B            <- fCrit*( 1 - r^2 )/resDF * r.factor
        bCI[1,1]     <- b*( sqrt(B+1) - sqrt(B) )
        bCI[1,2]     <- b*( sqrt(B+1) + sqrt(B) )
        varTest[1,1] <- vr[1,1] - 2*test.value*vr[1,2] + test.value^2*vr[2,2]
        varTest[1,2] <- vr[1,1] - test.value^2*vr[2,2]
        varTest[2,2] <- vr[1,1] + 2*test.value*vr[1,2] + test.value^2*vr[2,2]
    }
    else if ( (method==2) | (method=='MA') )
    {
        fac          <- vr[1,1] - vr[2,2]
        b            <- ( fac + sqrt( fac^2 + 4*vr[1,2]^2) )/2/vr[1,2]
        Q            <- fCrit*( vr[1,1]*vr[2,2] - vr[1,2]^2 )/resDF * r.factor
        bCI[1,1]     <- ( fac + sqrt( fac^2 + 4*vr[1,2]^2 - 4*Q) )/2/( vr[1,2] + sqrt(Q))
        bCI[1,2]     <- ( fac + sqrt( fac^2 + 4*vr[1,2]^2 - 4*Q) )/2/( vr[1,2] - sqrt(Q))
        if ( ( fac^2 + 4*vr[1,2]^2 - 4*Q) < 0 ) 
        {
            bCI[1,1] <- -Inf
            bCI[1,2] <-  Inf
        }
        varTest[1,1] <- vr[1,1] - 2*test.value*vr[1,2] + test.value^2*vr[2,2]
        varTest[1,2] <- -test.value^2*vr[1,2] + test.value*( vr[1,1] - vr[2,2] ) + vr[1,2]
        varTest[2,2] <- test.value^2*vr[1,1] + 2*test.value*vr[1,2] + vr[2,2]
    }
    else if ( (method==3) | (method=='lamest') )
    {   b=NA
        bCI[1,1:2]   <- NA
    }

     rTest  <- varTest[1,2] / sqrt( varTest[1,1] ) / sqrt( varTest[2,2] )
     F      <- rTest^2/(1 - rTest^2)/r.factor*(n-2)
     pValue <- 1 - pf( F, 1, resDF)

#      if ( is.null(data)==FALSE )
#      {
#         detach(data)
#      }

     list( F=F, r=rTest, p=pValue, test.value=test.value, b=b, ci=bCI )

}
