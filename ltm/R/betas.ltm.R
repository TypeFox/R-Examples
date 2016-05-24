betas.ltm <-
function (betas, constraint, p, q.) {
    if (!is.null(constraint)) {
        betas. <- matrix(0, p, q.)
        betas.[constraint[, 1:2, drop = FALSE]] <- constraint[, 3]
        betas.[-((constraint[, 2] - 1) * p + constraint[, 1])] <- betas
        betas.
    } else {
        matrix(betas, p)
    }
}
