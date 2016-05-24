hc2axes<- function (x) 
{
    A <- x$merge
    n <- nrow(A) + 1
    x.axis <- c()
    y.axis <- x$height
    x.tmp <- rep(0, 2)
    zz <- match(1:length(x$order), x$order)
    for (i in 1:(n - 1)) {
        ai <- A[i, 1]
        if (ai < 0) 
            x.tmp[1] <- zz[-ai]
        else x.tmp[1] <- x.axis[ai]
        ai <- A[i, 2]
        if (ai < 0) 
            x.tmp[2] <- zz[-ai]
        else x.tmp[2] <- x.axis[ai]
        x.axis[i] <- mean(x.tmp)
    }
    return(data.frame(x.axis = x.axis, y.axis = y.axis))
}
