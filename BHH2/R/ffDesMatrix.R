"ffDesMatrix" <-
function (k, gen = NULL) 
{
    N <- 2^k
    X <- matrix(NA, nrow = N, ncol = k)
    for (j in 1:k) X[, j] <- rep(sort(rep(c(-1, 1), N/2^j)), 
        2^(j - 1))
    X <- X[, ncol(X):1, drop=FALSE]
    if (is.null(gen)) 
        return(X)
    for (i in 1:length(gen)) {
        ind <- trunc(gen[[i]])
        if (any(abs(ind) > k)) 
            stop(paste("generator:", paste(ind[1], "=", paste(ind[-1], 
                collapse = "*")), "includes undefined columns"))
        x <- rep(sign(ind[1]), N)
        for (j in ind[-1]) x <- x * X[, j, drop=FALSE]
        X[, abs(ind[1])] <- x
    }
    X <- unique(X)
    X
}
