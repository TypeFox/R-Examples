LSM <-
function(x,p)
{
    x <- as.matrix(x)
    n <- nrow(x)
    xmat <- matrix(1,nrow=n-p,ncol=p+1)
    index <- p:(n-1)
    for(i in 1:p)
    {xmat[,i] <- x[index,1]
    index <- index-1
    }
   
    y <- x[(p+1):n,1]
    b <- solve( t(xmat) %*% xmat) %*% t(xmat) %*% y 
    #if (min(Mod(polyroot(c(1,b[1:p])))) <= 1) b[1:p] = ar.burg(y,aic=F,order.max=p)$ar
    e <- y - xmat %*% b
    cov1 <- sum(e^2)/(nrow(x)-length(b)) *  solve( t(xmat) %*% xmat)

return(list(coef=b,resid=e,covmat=cov1,xmat=xmat))
}
