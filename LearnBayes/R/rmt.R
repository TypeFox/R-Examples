rmt=function (n = 1, mean = rep(0, d), S, df = Inf) 
{
    d <- if (is.matrix(S)) 
        ncol(S)
    else 1
    if (df == Inf) 
        x <- 1
    else x <- rchisq(n, df)/df
    z <- rmnorm(n, rep(0, d), S)
    y <- t(mean + t(z/sqrt(x)))
    return(y)
}

