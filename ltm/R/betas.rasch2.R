betas.rasch2 <-
function (betas, constraint, p, IRT) {
    if (IRT) {
        if (!is.null(constraint)) {
            betas. <- numeric(p + 1)
            betas.[constraint[, 1]] <- constraint[, 2]
            betas.[-constraint[, 1]] <- betas
            cbind(- abs(betas.[p + 1]) * betas.[1:p], abs(betas.[p + 1]))
        } else {
            cbind(- abs(betas[p + 1]) * betas[1:p], abs(betas[p + 1]))
        }
    } else {
        if (!is.null(constraint)) {
            betas. <- numeric(p + 1)
            betas.[constraint[, 1]] <- constraint[, 2]
            betas.[-constraint[, 1]] <- betas
            cbind(betas.[1:p], abs(betas.[p + 1]))
        } else {
            cbind(betas[1:p], abs(betas[p + 1]))
        }
    }
}
