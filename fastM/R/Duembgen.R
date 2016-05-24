DUEMBGENshape <- function(X, nmax=500, eps = 1e-06, maxiter = 100, perm=FALSE)
    {
    Xnames <- colnames(X)

    X <- as.matrix(X)
    n <- nrow(X)
    if (perm) X <- X[sample(1:n),]

    if (n<=nmax)
        {
        res0 <- TYLERsymm1(X, eps = eps, maxiter = maxiter)
        Sigma <- res0$S
        rownames(Sigma) <- colnames(Sigma) <- Xnames
        } else {
        res0 <- TYLERsymm2(X, eps = eps, maxiter = maxiter)
        Sigma <- res0$S
        rownames(Sigma) <- colnames(Sigma) <- Xnames
        }
     
    if (res0$nG > eps) warning("DUEMBGENshape did not converge", call. = FALSE)

    return(list(Sigma=Sigma, iter=res0$iter))
    }



