"ImpulseCoefficientsARMA" <-
function(phi, theta, lag.max){
        p <- length(phi)
        q <- length(theta)
        if (p==0) return(c(-theta,rep(0,lag.max))[1:lag.max])
        r <- max(p,q)
        x <- numeric(lag.max + 1)
        t2<-numeric(r)
        t2[1:p]<-phi
        x[1] <- 1
        if (q>0) 
                x[2:(q+1)]<--theta
        for(i in 1:r) 
                x[i + 1] <- x[i + 1]+crossprod(t2[1:i], rev(x[1:i]))
        if(lag.max > r) 
                for(i in (r + 1):lag.max) 
                        x[i + 1] <- crossprod(phi, rev(x[(i-p+1):i]))            
        x
}

