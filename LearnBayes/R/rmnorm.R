rmnorm=function(n = 1, mean = rep(0, d), varcov) 
{
    d <- if (is.matrix(varcov)) 
        ncol(varcov)
    else 1
    z <- matrix(rnorm(n * d), n, d) %*% chol(varcov)
    y <- t(mean + t(z))
    return(y)
}
