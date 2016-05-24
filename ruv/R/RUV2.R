RUV2 <-
function (Y, X, ctl, k, Z = 1, eta = NULL, fullW = NULL, inputcheck = TRUE, 
    do_projectionplot = TRUE) 
{
    if (inputcheck) 
        inputcheck1(Y, X, Z, ctl)
    Y = RUV1(Y, eta, ctl)
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
    if (is.null(fullW)) {
        fullW = svd(Yc %*% t(Yc))$u[, 1:(m - p - q), drop = FALSE]
    }
    if (k > 0) {
        W = fullW[, 1:k, drop = FALSE]
        XZW = cbind(X, Z, W)
        if (do_projectionplot) {
            W0 = residop(W, X)
            bwx = solve(t(X) %*% X) %*% t(X) %*% W
            temp = svd(W0)
            vd = t((1/temp$d) * t(temp$v))
            W0 = W0 %*% vd
            bwx = bwx %*% vd
            W0Y = t(W0) %*% Y
            u = svd(W0Y %*% t(W0Y))$u
            W0 = W0 %*% u
            bwx = bwx %*% u
            projectionplotalpha = t(W0) %*% Y
            byx = solve(t(X) %*% X) %*% t(X) %*% Y
            projectionplotW = W0 + X %*% bwx
        }
        else {
            byx = bwx = projectionplotalpha = projectionplotW = NULL
        }
    }
    else {
        XZW = cbind(X, Z)
        W = alpha = byx = bwx = projectionplotW = projectionplotalpha = NULL
    }
    A = solve(t(XZW) %*% XZW)
    betagammaalphahat = A %*% t(XZW) %*% Y
    resids = Y - XZW %*% betagammaalphahat
    betahat = betagammaalphahat[1:p, , drop = FALSE]
    if (k > 0) 
        alpha = betagammaalphahat[(p + q + 1):(p + q + k), , 
            drop = FALSE]
    multiplier = as.matrix(diag(A)[1:p])
    df = m - p - q - k
    sigma2 = apply(resids^2, 2, sum)/df
    sigma2 = as.vector(sigma2)
    se = sqrt(multiplier %*% t(sigma2))
    tvals = betahat/se
    pvals = tvals
    for (i in 1:nrow(pvals)) pvals[i, ] = 2 * pt(-abs(tvals[i, 
        ]), df)
    return(list(betahat = betahat, sigma2 = sigma2, t = tvals, 
        p = pvals, multiplier = multiplier, df = df, W = W, alpha = alpha, 
        byx = byx, bwx = bwx, X = X, k = k, ctl = ctl, Z = Z, 
        fullW = fullW, projectionplotW = projectionplotW, projectionplotalpha = projectionplotalpha))
}
