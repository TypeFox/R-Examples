mAr.est <-
function (x, p, ...) 
{
    x = as.matrix(x)
    n = dim(x)[1]
    m = dim(x)[2]
    ne = n - p
    np = m * p + 1
    if (ne <= np) {
        stop("time series too short")
    }
    K = matrix(nrow = ne, ncol = np + m)
    K[, 1] = rep(1, ne)
    for (j in 1:p) {
        K[, seq(2 + m * (j - 1), 1 + m * j)] = data.matrix(x[seq(p - 
            j + 1, n - j), 1:m])
    }
    K[, seq(np + 1, np + m)] = data.matrix(x[seq(p + 1, n), 1:m])
    q = ncol(K)
    delta = (q^2 + q + 1) * (.Machine$double.eps)
    scale = sqrt(delta) * sqrt(apply(K^2, 2, sum))
    R = qr.R(qr((rbind(K, diag(scale)))), complete = TRUE)
    R22 = R[seq(np + 1, np + m), seq(np + 1, np + m)]
    logdp = 2 * log(abs(prod(diag(R22))))
    sbc = logdp/m - log(ne) * (ne - np)/ne
    R11 = R[seq(1, np), seq(1, np)]
    R12 = R[seq(1, np), seq(np + 1, np + m)]
    R22 = R[seq(np + 1, np + m), seq(np + 1, np + m)]
    R11[, 1] = R11[, 1] * max(scale[2:(np + m)])/scale[1]
    B = t(solve(R11) %*% R12)
    w = B[, 1] * max(scale[2:(np + m)])/scale[1]
    A = B[, 2:np]
    C = (t(R22) %*% R22)/(ne - np)
    t = seq(1, (n - p))
    res = matrix(nrow = (n - p), ncol = m)
    res[t, seq(1, m)] = x[t + p, ] - (rep(1, n - p) %*% t(w))
    for (j in seq(1, p)) {
        res[t, seq(1, m)] = res[t, seq(1, m)] - (x[(t - j + p), 
            ] %*% t(A[, seq(m * j - m + 1, j * m)]))
    }
    return(list(SBC = sbc, wHat = w, AHat = A, CHat = C, resid = res))
}
