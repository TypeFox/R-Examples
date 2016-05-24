thetas.tpm <-
function (thetas, type, constraint, p) {
    if (is.null(constraint)) {
        if (type == "rasch")
            matrix(c(thetas[-length(thetas)], abs(rep(thetas[length(thetas)], p))), p, 3)
        else
            matrix(thetas, p, 3)
    } else {
        thetas <- rep(thetas, length.out = 3 * p - nrow(constraint))
        thetas. <- matrix(0, p, 3)
        thetas.[constraint[, 1:2, drop = FALSE]] <- constraint[, 3]
        thetas.[-((constraint[, 2] - 1) * p + constraint[, 1])] <- thetas
        if (type == "rasch")
            thetas.[, 3] <- if (any(ind <- constraint[, 2] == 3)) abs(constraint[ind, 3]) else abs(thetas.[1, 3])
        thetas.
    }
}
