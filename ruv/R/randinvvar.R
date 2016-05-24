randinvvar <-
function (Y, ctl, XZ = NULL, eta = NULL, lambda = NULL, iterN = 1e+05) 
{
    getnormalX = function(m) {
        X = matrix(rnorm(m))
        X = X/sqrt(sum(X^2))
        return(X)
    }
    Y = RUV1(Y, eta, ctl)
    m = nrow(Y)
    n = ncol(Y)
    if (is.null(XZ)) {
        pq = 0
    }
    else if (length(XZ) == 1) {
        if (XZ == 1) {
            XZ = matrix(1, m, 1)
            pq = 1
        }
    }
    else {
        pq = ncol(XZ)
    }
    k = m - pq
    if (pq > 0) 
        temp = residop(Y, XZ)
    else temp = Y
    temp = svd(temp[, ctl] %*% t(temp[, ctl]))
    Y0 = t(temp$u[, 1:k, drop = FALSE]) %*% Y
    Y0cd = sqrt(temp$d[1:k])
    if (is.null(lambda)) 
        lambda = 0
    Y0cd2 = Y0cd^2 + lambda
    Ed = rep(0, k)
    for (i in 1:iterN) {
        Xstar = getnormalX(k)
        pseudoepsvect = Xstar/Y0cd2
        pseudoepsvect = pseudoepsvect/sqrt(sum(pseudoepsvect^2))
        Ed = Ed + pseudoepsvect^2
    }
    sigma2 = apply(Y0 * (as.vector(Ed) * Y0), 2, sum)/sum(Ed)
    df = sum(Ed)^2/sum(Ed^2)
    return(list(sigma2, df))
}
