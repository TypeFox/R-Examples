RUVinv <-
function (Y, X, ctl, Z = 1, eta = NULL, fullW0 = NULL, invsvd = NULL, 
    lambda = NULL, randomization = FALSE, iterN = 1e+05, inputcheck = TRUE) 
{
    if (inputcheck) 
        inputcheck1(Y, X, Z, ctl)
    Y = RUV1(Y, eta, ctl)
    if (is.null(lambda)) 
        do_ridge = FALSE
    else do_ridge = TRUE
    m = nrow(Y)
    n = ncol(Y)
    p = ncol(X)
    if (is.null(Z)) {
        q = 0
    }
    else if (length(Z) == 1) {
        if (Z == 1) {
            Z = matrix(1, m, 1)
            q = 1
        }
    }
    else {
        q = ncol(Z)
    }
    if (q > 0) {
        Y = residop(Y, Z)
        X = residop(X, Z)
    }
    Yc = Y[, ctl]
    Y0 = residop(Y, X)
    if (is.null(fullW0)) {
        fullW0 = svd(Y0 %*% t(Y0))$u[, 1:(m - p - q), drop = FALSE]
    }
    if (randomization) {
        if (do_ridge) 
            temp = randinvvar(Y, ctl, XZ = cbind(X, Z), lambda = lambda, 
                iterN = iterN)
        else temp = randinvvar(Y, ctl, XZ = cbind(X, Z), iterN = iterN)
    }
    else {
        if (do_ridge) 
            temp = invvar(Y, ctl, XZ = cbind(X, Z), lambda = lambda, 
                invsvd = invsvd)
        else temp = invvar(Y, ctl, XZ = cbind(X, Z), invsvd = invsvd)
    }
    sigma2 = temp[[1]]
    df = temp[[2]]
    sigma2 = as.vector(sigma2)
    k = m - p - q
    W0 = fullW0[, 1:k, drop = FALSE]
    alpha = solve(t(W0) %*% W0) %*% t(W0) %*% Y0
    byx = solve(t(X) %*% X) %*% t(X) %*% Y
    bycx = byx[, ctl, drop = FALSE]
    alphac = alpha[, ctl, drop = FALSE]
    if (do_ridge) 
        bwx = bycx %*% t(alphac) %*% solve(alphac %*% t(alphac) + 
            lambda * diag(nrow(alphac)))
    else bwx = bycx %*% t(alphac) %*% solve(alphac %*% t(alphac))
    W = W0 + X %*% bwx
    XZW = cbind(X, Z, W)
    A = solve(t(XZW) %*% XZW)
    betagammaalphahat = A %*% t(XZW) %*% Y
    betahat = betagammaalphahat[1:p, , drop = FALSE]
    multiplier = as.matrix(diag(A)[1:p])
    se = sqrt(multiplier %*% t(sigma2))
    tvals = betahat/se
    pvals = tvals
    for (i in 1:nrow(pvals)) pvals[i, ] = 2 * pt(-abs(tvals[i, 
        ]), df)
    return(list(betahat = betahat, sigma2 = sigma2, t = tvals, 
        p = pvals, multiplier = multiplier, df = df, W = W, alpha = alpha, 
        byx = byx, bwx = bwx, X = X, k = k, ctl = ctl, Z = Z, 
        fullW0 = fullW0, lambda = lambda, invsvd = temp$invsvd))
}
