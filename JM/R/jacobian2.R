jacobian2 <-
function (L, ncz) {
    ind <- which(lower.tri(matrix(0, ncz, ncz), TRUE), arr.ind = TRUE)
    dimnames(ind) <- NULL
    nind <- nrow(ind)
    id <- 1:nind
    rind <- which(ind[, 1] == ind[, 2])
    lind <- vector("list", length(rind))
    for (i in seq_along(rind)) {
        tt <- matrix(0, ncz - i + 1, ncz - i + 1)
        tt[lower.tri(tt, TRUE)] <- seq(rind[i], nind)
        tt <- tt + t(tt)
        diag(tt) <- diag(tt) / 2
        lind[[i]] <- tt
    }
    out <- matrix(0, nind, nind)
    for (g in 1:ncz) {
        gind <- id[g == ind[, 2]]
        vals <- L[gind]
        for (j in gind) {
            k <- which(j == gind)
            out[cbind(lind[[g]][k, ], j)] <- if (j %in% rind) vals[1] * vals else vals
        }
    }
    out[rind, ] <- 2 * out[rind, ]
    col.ind <- matrix(0, ncz, ncz)
    col.ind[lower.tri(col.ind, TRUE)] <- seq(1, length(L))
    col.ind <- t(col.ind)
    out[, col.ind[upper.tri(col.ind, TRUE)]]
}
