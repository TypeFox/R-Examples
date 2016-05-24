`nipals` <- function (Xtrain, Ytrain, Xtest = NULL, ncomp, tol = .Machine$double.eps^0.5, 
    weight = TRUE, beta = 0.1, ...) 
{
    Ytrain = as.matrix(Ytrain)
    n = dim(Xtrain)[1]
    p = dim(Xtrain)[2]
    m = dim(Ytrain)[2]
    R = P = matrix(, nrow = p, ncol = ncomp)
    Q = matrix(, nrow = ncomp, ncol = m)
    B = array(, dim = c(p, m, ncomp))
    T = U = matrix(, nrow = n, ncol = ncomp)
    tsqs = numeric(ncomp)
    fitted = residuals = array(, dim = c(n, m, ncomp))
    Ypred = array(, dim = c(dim(Xtest)[1], m, length(ncomp)))
    meanX = colMeans(Xtrain)
    X = Xtrain - rep(meanX, each = n)
    meanY = colMeans(Ytrain)
    Y = Ytrain - rep(meanY, each = n)
    if (weight == TRUE){
        wq = matrix(,n,1)
        for(i in 1:n){
            wq[i,] = beta*(1-beta)^(i-1)
        }
        weig = diag(rev(wq))
        X = weig%*%X
        Y = weig%*%as.matrix(Y)
    }
    Xtotvar = sum(Xtrain * Xtrain)
    for (a in 1:ncomp) {
        if (m == 1) {
            u = Y
        }
        else {
            u = Y[, which.max(colSums(Y * Y))]
            t.old = 0
        }
        repeat {
            r = crossprod(X, u)
            r = r/sqrt(c(crossprod(r)))
            t = X %*% r
            tsq = c(crossprod(t))
            t.t = t/tsq
            q = crossprod(Y, t.t)
            if (m == 1) 
                break
            if (sum(abs((t - t.old)/t), na.rm = T) < tol) 
                break
            else {
                u = Y %*% q/c(crossprod(q))
                t.old = t
            }
        }
        p = crossprod(X, t.t)
        X = X - t %*% t(p)
        Y = Y - t %*% t(q)
        R[, a] = r
        P[, a] = p
        Q[a, ] = q
        T[, a] = t
        U[, a] = u
        tsqs[a] = tsq
        fitted[, , a] = T[, 1:a] %*% Q[1:a, , drop = FALSE]
        residuals[, , a] = Y
    }
    if (ncomp == 1) {
        W = R
    }
    else {
        PR = crossprod(P, R)
        if (m == 1) {
            PRinv = diag(ncomp)
            bidiag = -PR[row(PR) == col(PR) - 1]
            for (a in 1:(ncomp - 1)) {
                PRinv[a, (a + 1):ncomp] = cumprod(bidiag[a:(ncomp - 
                  1)])
            }
        }
        else {
             PRinv = backsolve(PR, diag(ncomp))
        }
        W = R %*% PRinv
    }
    for (a in 1:ncomp) {
        B[, , a] = W[, 1:a, drop = FALSE] %*% Q[1:a, , drop = FALSE]
    }
    if (!is.null(Xtest)) {
        Xtest = scale(t(Xtest), scale = FALSE, center = meanX)
        Ypred = Xtest %*% B[, , length(ncomp)] + meanY
    }
    fitted = fitted + rep(meanY, each = n)
    class(T) = class(U) = "scores"
    class(P) = class(R) = class(Q) = "loadings"
    if (!is.null(Xtest)) {
        if (weight==TRUE){
            list(B = B, Ypred = Ypred, P = P, Q = t(Q), T = T, R = R, 
                 meanX = meanX, meanY = meanY, Yscores = U, projection = W, 
                 fitted.values = fitted, residuals = residuals, 
                 Xvar = colSums(P * P) * tsqs, Xtotvar = Xtotvar, weight = wq)
        }
        else
            list(B = B, Ypred = Ypred, P = P, Q = t(Q), T = T, R = R, 
                 meanX = meanX, meanY = meanY, Yscores = U, projection = W, 
                 fitted.values = fitted, residuals = residuals, 
                 Xvar = colSums(P * P) * tsqs, Xtotvar = Xtotvar)
    }
    else {
        if (weight==TRUE){
            list(B = B, P = P, Q = t(Q), T = T, R = R, meanX = meanX, 
                 meanY = meanY, Yscores = U, projection = W, fitted.values = fitted, 
                 residuals = residuals, Xvar = colSums(P * P) * tsqs, 
                 Xtotvar = Xtotvar, weight = wq)
        }
        else 
            list(B = B, P = P, Q = t(Q), T = T, R = R, meanX = meanX, 
                 meanY = meanY, Yscores = U, projection = W, fitted.values = fitted, 
                 residuals = residuals, Xvar = colSums(P * P) * tsqs, 
                 Xtotvar = Xtotvar)
    }
}
