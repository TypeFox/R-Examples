dvvtz<-function (v, z, dv, dz) 
{
    if (is.matrix(v) == FALSE) {
        v <- matrix(v, ncol = 1)
        dv <- array(dv, dim = c(1, nrow(dv), ncol(dv)))
    }
    k = ncol(v)
    p <- nrow(v)
    n<-dim(dv)[3]
    dummy <- matrix(0,dim(dv)[2],dim(dv)[3])
    for (i in 1:k) {
        D <- (v[, i] %*% t(z) + sum(v[, i] * z) * diag(p)) %*% 
            dv[i, , ] + v[, i] %*% t(v[, i]) %*% dz
        dummy <- dummy + D
    }
    return(dummy)
}
