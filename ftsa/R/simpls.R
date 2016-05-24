simpls = function(Xtrain, Ytrain, Xtest, ncomp = NULL, weight = FALSE, beta = 0.1)
{
    X = scale(Xtrain, center = TRUE, scale = FALSE)
    if (is.vector(Ytrain)){
       Ytrain = matrix(Ytrain, length(Ytrain), 1)
    }
    Y = scale(Ytrain, center = TRUE, scale = FALSE)
    meanX = apply(Xtrain, 2, mean)
    meanY = apply(Ytrain, 2, mean)
    n = dim(X)[1]
    p = dim(X)[2]
    m = dim(Y)[2]
    if (is.null(ncomp)){
       ncomp = min(n,p)
    }
    if (weight == TRUE){
       q = matrix(, n, 1)
       for(i in 1:n){
           q[i,] = beta * (1 - beta)^(i - 1)
       }
       weig = diag(rev(q))
       X = weig %*% X
       Y = weig %*% as.matrix(Y)
    }
    S = crossprod(X, Y)
    RR = matrix(0, ncol = max(ncomp), nrow = p)
    PP = matrix(0, ncol = max(ncomp), nrow = p)
    QQ = matrix(0, ncol = max(ncomp), nrow = m)
    TT = matrix(0, ncol = max(ncomp), nrow = n)
    VV = matrix(0, ncol = max(ncomp), nrow = p)
    UU = matrix(0, ncol = max(ncomp), nrow = n)
    B = array(0, c(dim(X)[2], dim(Y)[2], length(ncomp)))
    if (!is.null(Xtest)) {
        if (is.vector(Xtest)) {
            Xtest = matrix(Xtest, 1, length(Xtest))
        }
        Ypred = array(0, c(dim(Xtest)[1], m, length(ncomp)))
    }
    for (a in 1:max(ncomp)) {
        qq = svd(S)$v[, 1]
        rr = S %*% qq
        tt = scale(X %*% rr, scale = FALSE)
        tnorm = sqrt(sum(tt * tt))
        tt = tt / tnorm
        rr = rr / tnorm
        pp = crossprod(X, tt)
        qq = crossprod(Y, tt)
        uu = Y %*% qq
        vv = pp
        if (a > 1) {
            vv = vv - VV %*% crossprod(VV, pp)
            uu = uu - TT %*% crossprod(TT, uu)
        }
        vv = vv / sqrt(sum(vv * vv))
        S = S - vv %*% crossprod(vv, S)
        RR[, a] = rr
        TT[, a] = tt
        PP[, a] = pp
        QQ[, a] = qq
        VV[, a] = vv
        UU[, a] = uu
        if (!is.na(i <- match(a, ncomp))) {
            B[, , i] = RR[, 1:a, drop = FALSE] %*% t(QQ[, 1:a, 
                drop = FALSE])
            if (!is.null(Xtest)) {
                Xtest = scale(Xtest, scale = FALSE, center = meanX)
                Ypred[, , i] = Xtest %*% B[, , i] + meanY
            }
        }
    }
    if (length(ncomp) == 1) {
        B = B[, , 1]
    }
    if (!is.null(Xtest)) 
        if (weight==TRUE){
            list(B = B, Ypred = Ypred, P = PP, Q = QQ, T = TT, R = RR, weight = q)
        }
        else { 
            list(B = B, Ypred = Ypred, P = PP, Q = QQ, T = TT, R = RR)
        } 
    else 
        if (weight==TRUE){
            list(B = B, P = PP, Q = QQ, T = TT, R = RR, weight = q)
        }
        else {
            list(B = B, P = PP, Q = QQ, T = TT, R = RR)
        }
}
