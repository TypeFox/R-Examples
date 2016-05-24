deriv.D <-
function (D) {
    ncz <- nrow(D)
    ind <- which(lower.tri(D, TRUE), arr.ind = TRUE)
    dimnames(ind) <- NULL
    nind <- nrow(ind)
    svD <- solve(D)
    lapply(1:nind, function (x, ind) {
        mat <- matrix(0, ncz, ncz)
        ii <- ind[x, , drop = FALSE]
        mat[ii[1], ii[2]] <- mat[ii[2], ii[1]] <- 1
        mat
    },  ind = ind[, 2:1])
}
