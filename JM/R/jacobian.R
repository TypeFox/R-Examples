jacobian <-
function (theta) {
    k <- length(theta)
    etheta <- exp(theta)
    mat <- matrix(0, k, k)
    mat[, 1] <- rep(1, k)
    for (i in 2:k)
        mat[i:k, i] <- etheta[i]
    mat
}
