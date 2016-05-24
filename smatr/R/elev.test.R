
elev.test <- function( y, x, test.value=0, data=NULL, alpha=0.05, method="SMA", robust=FALSE, V=matrix(0,2,2) )
{
#     if ( is.null(data)==FALSE )
#     {
#         attach(data)
#     }

  
  if(!is.null(data))
    stop("'data' argument no longer supported.")
  
  
    iref <- ( is.na(x+y) == FALSE ) #to remove NA cases
    n    <- sum(iref)
    res.df <- n - 2
    fcrit  <- qf( 1-alpha, 1, res.df )
    dat    <- cbind( y[iref], x[iref] )
    if ( robust )
    {
		# get robust mean/var matrix:
		q     <- pchisq(3,2)
		S     <- huber.M(dat)
		means <- S$loc
		vr    <- ( S$cov - V) *(n-1)

		# get robust.factors (multipliers on variance matrix):
                rfac  <- robust.factor(dat,q)
 	        r.factor1 <- rfac[1]
                r.factor2 <- rfac[2]
    }
    else
    {
      r.factor1 <- 1
      r.factor2 <- 1 
      means    <- apply(dat,2,mean)
      vr <- ( var(dat) - V )*(n-1) 
    }	
    r      <- vr[1,2]/( ( vr[1,1]*vr[2,2] )^0.5 )

    if ( (method==0) | (method=="OLS") )
    {
        b       <- vr[1,2] / vr[2,2]
        var.res <- ( vr[1,1] - 2*b*vr[1,2] + b^2*vr[2,2] ) / res.df
        var.b   <- var.res / vr[2,2]
    }
    else if ( (method==1) | (method=="SMA") )
    {
        b       <- sign( vr[1,2] )*sqrt( vr[1,1] / vr[2,2] )
        var.res <- ( vr[1,1] - 2*b*vr[1,2] + b^2*vr[2,2] ) / res.df
        var.b   <- ( vr[1,1] - (vr[1,2]^2)/vr[2,2] ) / res.df / vr[2,2]
    }
    else if ( (method==2) | (method=="MA") )
    {
        fac     <- vr[1,1] - vr[2,2]
        b       <- ( fac + sqrt( fac^2 + 4*vr[1,2]^2) ) / 2 / vr[1,2]
        var.res <- ( vr[1,1] - 2*b*vr[1,2] + b^2*vr[2,2] ) / res.df
        var.fit <- ( b^2*vr[1,1] + 2*b*vr[1,2] + vr[2,2] ) / res.df
        var.b   <- 1 / ( var.res/var.fit + var.fit/var.res - 2)*( 1 + b^2 )^2 / res.df    # Use Fisher info
    }

    a        <- means[1] - b*means[2]
    var.a    <- var.res/n*r.factor2 + var.b*means[2]^2*r.factor1
    t        <- (a - test.value)/sqrt(var.a)
    pvalue   <- 2*pt( -abs(t), res.df )

#     if ( is.null(data)==FALSE )
#     {
#         detach(data)
#     }

    list( t=t, a=a, p=pvalue, a.ci=c( a-sqrt(var.a*fcrit), a+sqrt(var.a*fcrit) ), test.value=test.value )
}
