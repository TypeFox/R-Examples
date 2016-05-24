fullmodel.ssq <-
function (yX.data) 
{
    if (is.bma(yX.data)) {
        yX.data <- yX.data$X.data
    }
    y <- as.matrix(yX.data[, 1])
    X <- as.matrix(yX.data[, 2:ncol(yX.data)])
    N <- nrow(X)
    K = ncol(X)
    y.mean = mean(y)
    y <- y - matrix(y.mean, N, 1, byrow = TRUE)
    X.mean = colMeans(X)
    X <- X - matrix(X.mean, N, K, byrow = TRUE)
    Xqr <- qr(X)
    yty = as.numeric(crossprod(y))
    ymy = as.numeric(crossprod(qr.resid(Xqr, y)))
    ypy = as.numeric(crossprod(qr.fitted(Xqr, y)))
    R2 = ypy/yty
    return(list(R2 = R2, ymy = ymy, ypy = ypy, yty = yty, Fstat = (R2/(1 - 
        R2)) * (N - K - 1)/K))
}
