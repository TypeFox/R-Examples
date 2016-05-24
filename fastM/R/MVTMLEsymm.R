MVTMLEsymm <- function(X, nu = 1, nmax=500, eps = 1e-06, maxiter = 100, perm=FALSE)
    {
    Xnames <- colnames(X)

    X <- as.matrix(X)
    n <- nrow(X)
    if (perm) X <- X[sample(1:n),]

    if (n<=nmax)
        {
        res0 <- MVTMLEsymm1(X, nu = nu, eps = eps, maxiter = maxiter)
        Sigma <- res0$S
        rownames(Sigma) <- colnames(Sigma) <- Xnames
        } else {
        res0 <- MVTMLEsymm2(X, nu = nu, eps = eps, maxiter = maxiter)
        Sigma <- res0$S
        rownames(Sigma) <- colnames(Sigma) <- Xnames
        }
    if (res0$nG > eps) warning("MVTMLE did not converge", call. = FALSE)
    return(list(Sigma=Sigma, iter=res0$iter))
    }


