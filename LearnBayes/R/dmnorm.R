dmnorm=function (x, mean = rep(0, d), varcov, log = FALSE) 
{
    d <- if (is.matrix(varcov)) 
        ncol(varcov)
    else 1
    if (d > 1 & is.vector(x)) 
        x <- matrix(x, 1, d)
    n <- if (d == 1) 
        length(x)
    else nrow(x)
    X <- t(matrix(x, nrow = n, ncol = d)) - mean
    Q <- apply((solve(varcov) %*% X) * X, 2, sum)
    logDet <- sum(logb(abs(diag(qr(varcov)[[1]]))))
    logPDF <- as.vector(Q + d * logb(2 * pi) + logDet)/(-2)
    if (log) 
        logPDF
    else exp(logPDF)
}
