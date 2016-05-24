#lagIts-function----------------------------------------------------
lagIts <- function(x,k=1)
    {
    if (!inherits(x, "its")) stop("function is only valid for objects of class 'its'")
    lagmat <- x*NA
    dimnames(lagmat)[[2]] <- paste(dimnames(lagmat)[[2]],"lag",k)
    n <- dim(x)[1]
    if(k>0)
        {lagmat[(k+1):n,] <- x[1:(n-k),]}  else
        {lagmat[1:(n+k),] <- x[(1-k):n,]}
    return(lagmat)
    }
