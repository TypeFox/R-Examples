scio.cv <- function (X, lambda.max=1, alpha=0.95, cv.maxit=1e2, ...) {
    n <- nrow(X); p <- ncol(X)
    
    tr <- sample(1:n, round(n/2))
    S.tr <- cov(X[tr,])
    S.te <- cov(X[-tr,])
    lambda <- lambda.max
    lambda.cv <- lambda.max
    loss.min <- likelihood(S.te, scio(S.tr, lambda, ...)$w)
    for (itloss in 1:cv.maxit) {
        lambda <- lambda*alpha
        tmp <- likelihood(S.te, scio(S.tr, lambda, ...)$w)
        if (tmp <= loss.min) {
            loss.min <- tmp
            lambda.cv <- lambda
        } else {
            break
        }
    }

    if (itloss >= cv.maxit) warning("Maximum CV iterations exceeded! Consider increasing cv.maxit.")
    
    w <- scio(cov(X), lambda.cv, ...)$w
    return(list(w=w, lambda.cv=lambda.cv))
}

