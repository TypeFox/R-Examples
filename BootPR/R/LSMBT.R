LSMBT <-
function(x,p)
{
    x <- as.matrix(x)
    n <- nrow(x)
    xmat <- matrix(1,nrow=n-p,ncol=p+2)
    xmat[,p+2] <- 1:(n-p)
    index <- 2:(n-p+1)
    for(i in 1:p)
    {xmat[,i] <- x[index,1]
    index <- index+1
    }
   
    y <- x[1:(n-p),1]
    b <- solve( t(xmat) %*% xmat) %*% t(xmat) %*% y 
    e <- y - xmat %*% b

return(list(coef=b,resid=e))
}
