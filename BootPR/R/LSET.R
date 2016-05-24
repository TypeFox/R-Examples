LSET <-
function(x,p)
{
    x <- as.matrix(x)
    n <- nrow(x)
    xmat <- matrix(NA,nrow=n-p,ncol=p+2)
    xmat[,1] <- x[p:(n-1),1]
    if(p > 1)
    {
        z <- x[2:n,1] - x[1:(n-1),1]
        for(i in 1:(p-1))
            xmat[,i+1] <- z[(p-i):(n-1-i)]
    }

    xmat[,(p+1):(p+2)] <- cbind(rep(1,n-p),(p+1):n)
    
    y <- x[(p+1):n,1]
    b <- solve( t(xmat) %*% xmat) %*% t(xmat) %*% y 
    e <- y - xmat %*% b
    #return(list(coef=arlevel(b,p),resid=e))
    return(list(coef=b,resid=e))
}
