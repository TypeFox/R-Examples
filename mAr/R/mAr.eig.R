mAr.eig <-
function (A, C = NULL, ...) 
{
    m = dim(A)[1]
    p = (dim(A)[2])/m
    At = matrix(nrow = m * p, ncol = m * p)
    if (p == 1) 
        At = A
    else {
        At[seq(1, m), seq(1, m * p)] = A
        At[seq(m + 1, m * p), seq(1, m * p - m)] = diag(1, (p - 
            1) * m)
        At[seq(m + 1, m * p), seq(m * p - m + 1, m * p)] = 0
    }
    l = eigen(At)$values
    V = eigen(At)$vectors
    if ((any(Mod(l) > 1))) 
        warning("unstable AR model")
    a = matrix(nrow = 1, ncol = dim(V)[2])
    b = matrix(nrow = 1, ncol = dim(V)[2])
    St = matrix(nrow = dim(V)[2], ncol = dim(V)[2])
    for (j in seq(1, dim(V)[2])) {
        a = Re(V[, j])
        b = Im(V[, j])
        ph = 0.5 * atan(2 * sum(a * b)/(b %*% b - a %*% a))
        na = sqrt(sum((cos(ph) * a - sin(ph) * b)^2))
        nb = sqrt(sum((sin(ph) * a + sin(ph) * b)^2))
        if (nb > na && ph < 0) {
            ph = ph - pi/2
        }
        if (nb > na && ph > 0) {
            ph = ph + pi/2
        }
        St[, j] = V[, j] %*% exp((0 + (0+1i)) * ph)
    }
    S = St[seq(1 + (p - 1) * m, p * m), ]
    StInv = solve(St)[, seq(1, m)]
    tau = matrix(nrow = 1, ncol = m * p)
    per = matrix(nrow = 1, ncol = m * p)
    exctn = matrix(nrow = 1, ncol = m * p)
    for (i in seq(1, m * p)) {
        tau[i] = -2/log((abs(l[i]))^2)
        a = Re(l[i])
        b = Im(l[i])
        if (identical(b, 0) && a >= 0) {
            per[i] = Inf
        }
        if (identical(b, 0) && a < 0) {
            per[i] = 2
        }
        else {
            per[i] = 2 * pi/abs(atan2(b, a))
        }
    }
    M = cbind(periods = as.vector(per), dampTime = as.vector(tau))
    result = M[c(order(M[, 2], decreasing = TRUE)), ]
    return(list(modes = result, eigv = S))
}
