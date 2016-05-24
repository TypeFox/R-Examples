
b.com.est <- function( z, n, method, lambda=1, res.df)
{
    zcom    <- t(z) %*% ( n - 1 )
    zlam    <- zcom
    iter    <- 100
    bi      <- matrix( 0, iter, 1 )
    bchange <- 1
    i       <- 1

    while ( abs(bchange) > 10^(-6) )
    {
        i  <- i+1;
        if ( i > 100 ) 
        {
            warning("no convergence!")
		i <- i-1 #so that i=100
            break
        }

        if ( (method==1) | (method=="SMA") )
        {
            bi[i] <- sign( zcom[3] )*sqrt( zcom[1]/zcom[2] )
            l1    <- ( bi[i]^2*z[,2] + 2*bi[i]*z[,3] + z[,1] )/2/abs(bi[i])
            l2    <- ( bi[i]^2*z[,2] - 2*bi[i]*z[,3] + z[,1] )/2/abs(bi[i])
            wts   <- n * (1/l2 + 1/l1)
            wts[n==1] <- 0 #to avoid errors with n==1
            zcom  <- t(z) %*% wts
        }
        else if ( (method==2) | (method=="MA") )
        {
            fac   <- zcom[1] - lambda*zcom[2]
            bi[i] <- ( fac + sqrt( fac^2 + 4*lambda*zcom[3]^2 ) ) /2/zcom[3]
            l1    <- ( lambda^2*z[,2] + 2*lambda*bi[i]*z[,3] + bi[i]^2*z[,1] ) / (lambda+bi[i]^2)
            l2    <- ( bi[i]^2*z[,2] - 2*bi[i]*z[,3] + z[,1] )/ ( lambda + bi[i]^2 )
            wts   <- n * (1/l2 - lambda/l1)
            wts[n==1] <- 0 #to avoid errors with n==1
            zcom  <- t(z) %*% wts
        }
        else if ( (method==3) | (method=="lamest") )
        {
            fac   <- zcom[1] - lambda*zcom[2]
    	      bi[i] <- ( fac + sqrt( fac^2 + 4*lambda*zcom[3]^2 ) ) /2/zcom[3]
            lambda<- abs( ( bi[i]^2*zlam[3] - bi[i]*zlam[1] ) / (zlam[3] - bi[i]*zlam[2]))
            l1    <- ( lambda^2*z[,2] + 2*lambda*bi[i]*z[,3] + bi[i]^2*z[,1] ) / (lambda+bi[i]^2)
            l2    <- ( bi[i]^2*z[,2] - 2*bi[i]*z[,3] + z[,1] )/ ( lambda + bi[i]^2 )
            wts   <- n * (1/l2 - lambda/l1)
            wts[n==1] <- 0 #to avoid errors with n==1
            zcom  <- t(z) %*% wts
            wtLam <- n / l1
            wtLam[n==1] <- 0 #to avoid errors with n==1
            zlam  <- t(z) %*% wtLam
        }
        else 
        {
            stop("No such method")
        }

        bchange <- bi[i] - bi[i-1]
        b       <- bi[i]
    }
    if ( (method==1) | (method=='SMA') ) { lambda <- b^2 }

    bi=bi[2:i] #ignoring first entry, which was 0.

    list( b=b, bi=bi, l1=l1, l2=l2, lambda=lambda )
}



