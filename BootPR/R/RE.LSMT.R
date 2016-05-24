RE.LSMT <-
function(x,p,b)
{
    x <- as.matrix(x)
    n <- nrow(x)
    xmat <- matrix(1,nrow=n-p,ncol=p+2)
    xmat[,p+2] <- (p+1):n
    index <- p:(n-1)
    for(i in 1:p)
    {xmat[,i] <- x[index,1]
    index <- index-1
    }
   
    y <- x[(p+1):n,1]
    x1 <- as.matrix(xmat[,1:p])
    x2 <- xmat[,(p+1):(p+2)]
    tem1 <- y-x1 %*% b[1:p,1]
    tem2 <- solve( t(x2) %*% x2) %*% t(x2) %*% tem1 
 return(tem2)
 }
