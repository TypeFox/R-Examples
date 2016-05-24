RESID <-
function(x,b)
{
    x <- as.matrix(x)
    n <- nrow(x)
    p <- length(b)-1
    xmat <- matrix(1,nrow=n-p,ncol=p+1)
    index <- p:(n-1)
    for(i in 1:p)
    {xmat[,i] <- x[index,1]
    index <- index-1
    }
   
    y <- x[(p+1):n,1]
    e <- y - xmat %*% b
    return(e-mean(e))
}
