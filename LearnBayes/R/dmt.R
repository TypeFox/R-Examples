dmt=function (x, mean = rep(0, d), S, df = Inf, log = FALSE) 
{
    if (df == Inf) 
        return(dmnorm(x, mean, S, log = log))
    d <- if (is.matrix(S)) 
        ncol(S)
    else 1
    if (d > 1 & is.vector(x)) 
        x <- matrix(x, 1, d)
    n <- if (d == 1) 
        length(x)
    else nrow(x)
    X <- t(matrix(x, nrow = n, ncol = d)) - mean
    Q <- apply((solve(S) %*% X) * X, 2, sum)
    logDet <- sum(logb(abs(diag(qr(S)$qr))))
    logPDF <- (lgamma((df + d)/2) - 0.5 * (d * logb(pi * df) + 
        logDet) - lgamma(df/2) - 0.5 * (df + d) * logb(1 + Q/df))
    if (log) 
        logPDF
    else exp(logPDF)
}

