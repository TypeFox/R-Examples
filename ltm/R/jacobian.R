jacobian <-
function (thetas, constrained, ind1, ind2, p) {
    lis.mat <- vector("list", p)
    for (j in 1:p) {
        theta <- if (constrained) thetas[seq(ind1[j], ind2[j])] else thetas[seq(ind1[j], ind2[j] - 1)]
        nth <- length(theta)
        etheta <- exp(theta)
        mat <- matrix(0, nth, nth)
        mat[, 1] <- rep(1, nth)
        if (nth > 1) {
            for (i in 2:nth)
                mat[i:nth, i] <- etheta[i]
        }
        lis.mat[[j]] <- mat
    }
    lis.mat
}
