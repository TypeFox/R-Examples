mAr.pca <-
function (x, p, k = dim(x)[2], ...) 
{
    require(mAr)
    n = dim(x)[1]
    m = dim(x)[2]
    if (k <= 1) {
        stop("EOF subspace dimension <=1")
    }
    xcentred = matrix(nrow = n, ncol = m)
    xscaled = matrix(nrow = n, ncol = m)
    for (i in 1:m) {
        xcentred[, i] = x[, i] - apply(x, 2, mean)[i]
    }
    for (i in 1:m) {
        xscaled[, i] = xcentred[, i]/apply(xcentred, 2, sd)[i]
    }
    d = svd(xscaled)$d
    fraction.variance = cumsum(d^2)/sum(d^2)
    D = diag(svd(xscaled)$d[1:k])
    U = svd(xscaled)$u[, 1:k]
    V = svd(xscaled)$v[, 1:k]
    scores = U %*% D
    loadings = V
    mArModel = mAr.est(scores, p)
    A = mArModel$AHat
    C = mArModel$CHat
    SBC = mArModel$SBC
    modes = mAr.eig(A)$modes
    P = V %*% mAr.eig(A, C)$eigv
    return(list(order = p, SBC = SBC, fraction.variance = fraction.variance[1:k], 
        resid = mArModel$resid, eigv = P, modes = modes, scores = scores))
}
